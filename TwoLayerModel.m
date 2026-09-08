function [temperature,displacement] = TwoLayerModel(t,Lambda,p,S)
% Calculate film-surface temperature and displacement responses
% Inputs:
%   t - Pump-relative time samples [s]
%   Lambda - Grating spacing [m]
%   p - [alpha_f, alpha_s, R]
%   S - Known material and Fourier-reconstruction settings
% Outputs:
%   temperature - Surface temperature response at each time
%   displacement - Surface-normal displacement response at each time

%% Validate and orient model inputs
    t = t(:);
    p = p(:).';
    if ~isscalar(Lambda) || ~isfinite(Lambda) || Lambda <= 0
        error("TGS:Grating","Lambda must be a positive finite scalar.");
    end
    if numel(p) ~= 3 || any(~isfinite(p)) || any(p <= 0)
        error("TGS:Parameters", "p must contain positive finite [alpha_f, alpha_s, R].");
    end
    if any(~isfinite(t))
        error("TGS:Time","Time samples must be finite.");
    end
    if ~isscalar(S.Nw) || S.Nw < 4 || mod(S.Nw,2) ~= 0
        error("TGS:FourierGrid","S.Nw must be an even integer of at least four.");
    end

%% Calculate material and coupling constants
    alpha_f = p(1);
    alpha_s = p(2);
    R = p(3);
    L = S.L;
    q = 2*pi/Lambda;
    filmConductivity = S.Cf*alpha_f;
    substrateConductivity = S.Cs*alpha_s;
    Lame_f = S.Ef/(2*(1+S.nuf));
    Lame_s = S.Es/(2*(1+S.nus));
    gamma_f = S.alphathf*(1+S.nuf)/(1-S.nuf);
    gamma_s = S.alphaths*(1+S.nus)/(1-S.nus);

%% Build the symmetric angular-frequency grid
    positiveCount = S.Nw/2;
    omegaPositive = logspace(log10(S.wmin),log10(S.wmax),positiveCount);
    omega = [-fliplr(omegaPositive),omegaPositive];
    temperaturePositive = complex(zeros(size(omegaPositive)));
    displacementPositive = complex(zeros(size(omegaPositive)));

%% Assemble the elastic boundary-condition system
    exp2qL = exp(2*q*L);
    shearRatio = Lame_s/Lame_f;
    elasticMatrix = [ 1, 2*S.nuf-2,      1,      2-2*S.nuf,               0,           0; ...
                      1, 2*S.nuf-1,     -1,      2*S.nuf-1,               0,           0; ...
                      1, q*L-3+4*S.nuf, -exp2qL, exp2qL*(4*S.nuf-q*L-3), -1,           3-4*S.nus-q*L; ...
                      1, q*L,            exp2qL, q*L*exp2qL,             -1,          -q*L; ...
                      1, q*L-2+2*S.nuf,  exp2qL, exp2qL*(q*L+2-2*S.nuf), -shearRatio,  (2-2*S.nus-q*L)*shearRatio; ...
                      1, q*L-1+2*S.nuf, -exp2qL, exp2qL*(2*S.nuf-q*L-1), -shearRatio,  (1-2*S.nus-q*L)*shearRatio];
    elasticSolver = decomposition(elasticMatrix,"lu"); % factors matrix using LU decomposition

%% Solve the coupled response at each positive frequency
    for frequencyIndex = 1:positiveCount
        currentOmega = omegaPositive(frequencyIndex);
        beta_f = sqrt(q^2+1i*currentOmega/alpha_f);
        beta_s = sqrt(q^2+1i*currentOmega/alpha_s);
        r_f = 1i*currentOmega/alpha_f;
        r_s = 1i*currentOmega/alpha_s;

        conductivityRatio = substrateConductivity*beta_s/(filmConductivity*beta_f);
        interfaceFactor = 1 - conductivityRatio + substrateConductivity*beta_s*R;
        thermalDenominator = 2*conductivityRatio + interfaceFactor*(1-exp(-2*beta_f*L));
        commonFactor = S.Q0/(filmConductivity*beta_f*thermalDenominator);
        interfaceDecay = exp((q-beta_f)*L);
        A_f = commonFactor*(2*conductivityRatio+interfaceFactor);
        B_f = commonFactor*interfaceFactor*exp(-2*beta_f*L);
        A_f_interface = commonFactor* (2*conductivityRatio+interfaceFactor)*interfaceDecay;
        B_f_interface = commonFactor*interfaceFactor*interfaceDecay;
        A_s_interface = 2*commonFactor*interfaceDecay;

        elasticForcing = [ gamma_f*beta_f*(A_f-B_f)/r_f; ...
                           gamma_f*q*(A_f+B_f)/r_f; ...
                           gamma_f*q*(A_f_interface+B_f_interface)/r_f-gamma_s*q*A_s_interface/r_s; ...
                           gamma_f*beta_f*(A_f_interface-B_f_interface)/r_f-gamma_s*beta_s*A_s_interface/r_s; ...
                           gamma_f*beta_f*(A_f_interface-B_f_interface)/r_f-shearRatio*gamma_s*beta_s*A_s_interface/r_s; ...
                           gamma_f*q*(A_f_interface+B_f_interface)/r_f-shearRatio*gamma_s*q*A_s_interface/r_s];
        elasticCoefficients = elasticSolver\elasticForcing;

        temperaturePositive(frequencyIndex) = A_f + B_f;
        displacementPositive(frequencyIndex) = gamma_f*beta_f*(B_f-A_f)/r_f + elasticCoefficients(1)+elasticCoefficients(3);
    end

%% Reconstruct the time-domain responses
    temperatureSpectrum = [conj(fliplr(temperaturePositive)), temperaturePositive];
    displacementSpectrum = [conj(fliplr(displacementPositive)), displacementPositive];
    temperature = real(inverseFourierPiecewiseLinear( omega,temperatureSpectrum,t))/(2*pi);
    displacement = real(inverseFourierPiecewiseLinear( omega,displacementSpectrum,t))/(2*pi);
end

function timeSignal = inverseFourierPiecewiseLinear(omega,spectrum,t)
% Integrate a linearly interpolated spectrum

    omega = omega(:).';
    spectrum = spectrum(:).';
    t = t(:);
    deltaOmega = diff(omega);
    spectralSlope = diff(spectrum)./deltaOmega;
    timeSignal = complex(zeros(size(t)));

    for timeIndex = 1:numel(t)
        currentTime = t(timeIndex);
        scaledInterval = currentTime*deltaOmega;
        useSeries = abs(scaledInterval) < 1e-4;
        integralConstant = complex(zeros(size(deltaOmega)));
        integralLinear = complex(zeros(size(deltaOmega)));

        % Use closed forms away from zero and series expansions near zero.
        interval = deltaOmega(~useSeries);
        phase = exp(1i*currentTime*interval);
        integralConstant(~useSeries) = (phase-1)/(1i*currentTime);
        integralLinear(~useSeries) = phase.*(interval/(1i*currentTime)+1/currentTime^2)- 1/currentTime^2;

        interval = deltaOmega(useSeries);
        integralConstant(useSeries) = interval + 1i*currentTime*interval.^2/2-currentTime^2*interval.^3/6;
        integralLinear(useSeries) = interval.^2/2 + 1i*currentTime*interval.^3/3-currentTime^2*interval.^4/8;

        timeSignal(timeIndex) = sum( exp(1i*currentTime*omega(1:end-1)).* (spectrum(1:end-1).* integralConstant + spectralSlope.*integralLinear));
    end
end
