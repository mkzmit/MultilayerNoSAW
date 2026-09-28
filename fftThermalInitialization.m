function [initializationTraces,diagnostics] = fftThermalInitialization(traces,S)
% Build a thermal-only initialization trace by FFT/IFFT band removal.
% No acoustic waveform is fitted. The measured SAW band is detected from a
% thermal-background residual, reconstructed with a conjugate-symmetric FFT
% mask, and subtracted only in the returned initialization trace.

    if isempty(traces) || ~isstruct(traces)
        error("TGS:SAWTrace", ...
            "At least one prepared trace structure is required.");
    end

    if isfield(traces,"fftInitialized")
        alreadyInitialized = arrayfun(@(trace) ...
            isequal(trace.fftInitialized,true),traces);
        if all(alreadyInitialized)
            if ~isfield(traces,"fftInitialization")
                error("TGS:SAWTrace", ...
                    "An FFT-initialized trace must retain its diagnostics.");
            end
            initializationTraces = traces;
            diagnostics = reshape([traces.fftInitialization],[],1);
            return
        elseif any(alreadyInitialized)
            error("TGS:SAWTrace", ...
                "Do not mix raw and already FFT-initialized traces in one call.");
        end
    end

    traceCount = numel(traces);
    sourceTime = cell(traceCount,1);
    sourceSignal = cell(traceCount,1);
    sourceResidual = cell(traceCount,1);
    backgroundSource = strings(traceCount,1);
    sampleInterval = zeros(traceCount,1);
    sampleRate = zeros(traceCount,1);
    recordDuration = zeros(traceCount,1);
    lambda = zeros(traceCount,1);

    for traceIndex = 1:traceCount
        [sourceTime{traceIndex},sourceSignal{traceIndex}] = ...
            fullResolutionTrace(traces(traceIndex));
        sampleInterval(traceIndex) = validateFFTGrid( ...
            sourceTime{traceIndex},sourceSignal{traceIndex},S);
        sampleRate(traceIndex) = 1/sampleInterval(traceIndex);
        recordDuration(traceIndex) = ...
            numel(sourceTime{traceIndex})*sampleInterval(traceIndex);

        if ~isfield(traces(traceIndex),"Lambda") || ...
                ~isscalar(traces(traceIndex).Lambda) || ...
                ~isfinite(traces(traceIndex).Lambda) || ...
                traces(traceIndex).Lambda <= 0
            error("TGS:SAWWavelength", ...
                "Every trace must contain a positive finite Lambda in meters.");
        end
        lambda(traceIndex) = double(traces(traceIndex).Lambda);
        [sourceResidual{traceIndex},backgroundSource(traceIndex)] = ...
            thermalOnlyResidual(sourceTime{traceIndex}, ...
            sourceSignal{traceIndex},lambda(traceIndex),S);
    end

    settings = sawSettings(S,lambda,sampleRate,recordDuration);
    if isfinite(settings.configuredFrequencyHz)
        sawFrequencyHz = settings.configuredFrequencyHz;
        peakRatio = Inf;
        requiredPeakRatio = NaN;
        detected = true;
        detectionSource = "configured";
    else
        [sawFrequencyHz,peakRatio,requiredPeakRatio,detected] = ...
            detectSAWFrequency( ...
            sourceTime,sourceResidual,sampleInterval,settings);
        detectionSource = "fft_derivative_psd";
    end

    if detected
        halfBandwidthHz = settings.halfBandwidthHz;
        if ~isfinite(halfBandwidthHz)
            halfBandwidthHz = max(settings.bandFraction*sawFrequencyHz, ...
                settings.minimumHalfBandwidthHz);
        end
    else
        sawFrequencyHz = NaN;
        halfBandwidthHz = NaN;
        warning("TGS:SAWNotDetected", ...
            "No FFT peak in %.4g-%.4g MHz exceeded the required " + ...
            "peak ratio of %.3g; the initialization traces were left unchanged.", ...
            settings.searchBandHz(1)/1e6,settings.searchBandHz(2)/1e6, ...
            requiredPeakRatio);
    end

    if detected && halfBandwidthHz < settings.nativeResolutionHz
        error("TGS:SAWBandwidth", ...
            "The SAW half-bandwidth must span at least one native FFT bin.");
    end
    if detected && (halfBandwidthHz >= sawFrequencyHz || ...
            sawFrequencyHz+halfBandwidthHz >= min(sampleRate)/2)
        error("TGS:SAWBandwidth", ...
            "The SAW band must remain strictly between zero and Nyquist.");
    end

    diagnosticTemplate = struct( ...
        "method","fft_band_ifft", ...
        "detected",false, ...
        "frequencyHz",NaN, ...
        "searchBandHz",settings.searchBandHz, ...
        "halfBandwidthHz",NaN, ...
        "peakRatio",peakRatio, ...
        "requiredPeakRatio",requiredPeakRatio, ...
        "detectionSource",detectionSource, ...
        "backgroundSource","", ...
        "bandPowerBefore",NaN, ...
        "bandPowerAfter",NaN, ...
        "bandPowerRemovedFraction",NaN, ...
        "sawRMS",0);
    diagnostics = repmat(diagnosticTemplate,traceCount,1);
    initializationTraces = traces;

    for traceIndex = 1:traceCount
        if detected
            [sawTermFull,powerBefore,powerAfter] = extractSAWBand( ...
                sourceTime{traceIndex},sourceResidual{traceIndex}, ...
                sampleInterval(traceIndex),sawFrequencyHz, ...
                halfBandwidthHz,settings);
        else
            sawTermFull = zeros(size(sourceSignal{traceIndex}));
            powerBefore = NaN;
            powerAfter = NaN;
        end

        cleanedFull = sourceSignal{traceIndex}-sawTermFull;
        fitTime = double(traces(traceIndex).t(:));
        fitRaw = double(traces(traceIndex).y(:));
        sawTerm = interp1(sourceTime{traceIndex},sawTermFull, ...
            fitTime,"linear");
        if any(~isfinite(sawTerm))
            error("TGS:SAWInterpolation", ...
                "The fit samples must lie inside the full FFT time grid.");
        end

        diagnostics(traceIndex).detected = detected;
        diagnostics(traceIndex).frequencyHz = sawFrequencyHz;
        diagnostics(traceIndex).halfBandwidthHz = halfBandwidthHz;
        diagnostics(traceIndex).backgroundSource = backgroundSource(traceIndex);
        diagnostics(traceIndex).bandPowerBefore = powerBefore;
        diagnostics(traceIndex).bandPowerAfter = powerAfter;
        if isfinite(powerBefore) && powerBefore > 0
            diagnostics(traceIndex).bandPowerRemovedFraction = ...
                1-powerAfter/powerBefore;
        end
        postPump = sourceTime{traceIndex} >= 0;
        if any(postPump)
            diagnostics(traceIndex).sawRMS = ...
                sqrt(mean(abs(sawTermFull(postPump)).^2));
        end

        initializationTraces(traceIndex).t = fitTime;
        initializationTraces(traceIndex).yRaw = fitRaw;
        initializationTraces(traceIndex).y = fitRaw-sawTerm;
        initializationTraces(traceIndex).sawTerm = sawTerm;
        initializationTraces(traceIndex).tFull = sourceTime{traceIndex};
        initializationTraces(traceIndex).yFullRaw = sourceSignal{traceIndex};
        initializationTraces(traceIndex).yFull = cleanedFull;
        initializationTraces(traceIndex).sawTermFull = sawTermFull;
        initializationTraces(traceIndex).fftInitialization = diagnostics(traceIndex);
        initializationTraces(traceIndex).fftInitialized = true;
    end
end

function [t,y] = fullResolutionTrace(trace)
% Prefer the uniform acquisition grid preserved by prepareTGSData.

    hasFullTime = isfield(trace,"tFull") && ~isempty(trace.tFull);
    hasFullSignal = isfield(trace,"yFull") && ~isempty(trace.yFull);
    if hasFullTime ~= hasFullSignal
        error("TGS:SAWTrace", ...
            "tFull and yFull must either both be present or both be absent.");
    end

    if hasFullTime
        t = double(trace.tFull(:));
        y = double(trace.yFull(:));
    else
        if ~isfield(trace,"t") || ~isfield(trace,"y")
            error("TGS:SAWTrace", "Each trace must contain t and y arrays.");
        end
        t = double(trace.t(:));
        y = double(trace.y(:));
    end
end

function dt = validateFFTGrid(t,y,S)
% Require a finite, increasing, nearly uniform time grid for the FFT.

    if numel(t) ~= numel(y) || numel(t) < 16 || ...
            any(~isfinite(t)) || any(~isfinite(y))
        error("TGS:FFTGrid", ...
            "FFT thermal initialization requires at least 16 paired finite samples.");
    end

    timeStep = diff(t);
    if any(timeStep <= 0)
        error("TGS:FFTGrid", ...
            "FFT thermal initialization requires strictly increasing time samples.");
    end

    uniformTolerance = 1e-5;
    if isfield(S,"sawUniformTolerance")
        uniformTolerance = double(S.sawUniformTolerance);
    end
    if ~isscalar(uniformTolerance) || ~isfinite(uniformTolerance) || ...
            uniformTolerance < 0
        error("TGS:SAWSettings", ...
            "S.sawUniformTolerance must be a nonnegative finite scalar.");
    end

    dt = median(timeStep);
    relativeJitter = max(abs(timeStep-dt))/dt;
    if relativeJitter > uniformTolerance
        error("TGS:FFTGrid", ...
            "FFT thermal initialization requires a uniform time grid " + ...
            "(relative step variation %.3g exceeds %.3g). " + ...
            "Provide uniform tFull/yFull arrays.", ...
            relativeJitter,uniformTolerance);
    end
end

function [residual,source] = thermalOnlyResidual(t,y,lambda,S)
% Remove a coarse thermal-only background before the acoustic FFT.
% This follows the reference workflow's background subtraction without
% adding any acoustic component or acoustic parameter to the model.

    required = ["p0","L","Cf","Cs","Ef","Es","nuf","nus", ...
        "alphathf","alphaths","Nw","wmin","wmax","Q0"];
    hasThermalModel = all(isfield(S,required));
    postPump = t >= 0;

    if hasThermalModel && nnz(postPump) >= 3
        [temperature,displacement] = TwoLayerModel( ...
            t(postPump),lambda,S.p0,S);
        design = [temperature(:),displacement(:),ones(nnz(postPump),1)];
        columnScale = vecnorm(design,2,1);
        if any(~isfinite(columnScale)) || any(columnScale == 0)
            error("TGS:SAWBackground", ...
                "The coarse thermal-only background is singular.");
        end
        scaledDesign = design./columnScale;
        scaledCoefficient = scaledDesign\y(postPump);
        coefficient = scaledCoefficient./columnScale.';

        background = zeros(size(y));
        background(postPump) = design*coefficient;
        if any(~postPump)
            background(~postPump) = mean(y(~postPump));
        end
        residual = y-background;
        source = "thermal_model";
    else
        residual = y;
        source = "endpoint_detrend";
    end

    residual = endpointDetrend(residual);
end

function settings = sawSettings(S,lambda,sampleRate,recordDuration)
% Validate settings and construct a shared physically constrained search band.

    nativeResolution = max(1./recordDuration);
    positiveLimit = min(sampleRate/2-1./recordDuration);

    configuredFrequencyHz = NaN;
    if isfield(S,"sawFrequencyHz") && ~isempty(S.sawFrequencyHz)
        configuredFrequencyHz = double(S.sawFrequencyHz);
        if ~isscalar(configuredFrequencyHz) || ...
                ~isfinite(configuredFrequencyHz) || ...
                configuredFrequencyHz <= 0 || ...
                configuredFrequencyHz >= positiveLimit
            error("TGS:SAWSettings", ...
                "S.sawFrequencyHz must be positive and below Nyquist.");
        end
    end

    if isfield(S,"sawFrequencyRange")
        frequencyRange = double(S.sawFrequencyRange(:).');
        if numel(frequencyRange) ~= 2 || ...
                any(~isfinite(frequencyRange)) || ...
                frequencyRange(1) <= 0 || ...
                frequencyRange(1) >= frequencyRange(2)
            error("TGS:SAWSettings", ...
                "S.sawFrequencyRange must be [minimum, maximum] in Hz.");
        end
        searchBandHz = frequencyRange;
    elseif isfield(S,"sawVelocityRange")
        velocityRange = double(S.sawVelocityRange(:).');
        if numel(velocityRange) ~= 2 || ...
                any(~isfinite(velocityRange)) || ...
                velocityRange(1) <= 0 || ...
                velocityRange(1) >= velocityRange(2)
            error("TGS:SAWSettings", ...
                "S.sawVelocityRange must be [minimum, maximum] in m/s.");
        end
        perTraceBand = velocityRange./lambda;
        searchBandHz = [max(perTraceBand(:,1)),min(perTraceBand(:,2))];
    elseif isfinite(configuredFrequencyHz)
        searchBandHz = [configuredFrequencyHz,configuredFrequencyHz];
    else
        error("TGS:SAWSettings", ...
            "Specify S.sawVelocityRange in m/s or " + ...
            "S.sawFrequencyRange in Hz for FFT peak detection.");
    end

    if ~isfinite(configuredFrequencyHz)
        searchBandHz(2) = min(searchBandHz(2),positiveLimit);
        if searchBandHz(1) >= searchBandHz(2) || ...
                diff(searchBandHz) < 4*nativeResolution
            error("TGS:SAWSearchBand", ...
                "The SAW search band does not contain enough resolvable " + ...
                "frequencies below the Nyquist limit (%.4g MHz).", ...
                positiveLimit/1e6);
        end
    end

    settings.searchBandHz = searchBandHz;
    settings.zeroPadFactor = 8;
    if isfield(S,"sawZeroPadFactor")
        settings.zeroPadFactor = double(S.sawZeroPadFactor);
    end
    if ~isscalar(settings.zeroPadFactor) || ...
            ~isfinite(settings.zeroPadFactor) || ...
            settings.zeroPadFactor < 1 || ...
            settings.zeroPadFactor ~= fix(settings.zeroPadFactor)
        error("TGS:SAWSettings", ...
            "S.sawZeroPadFactor must be a positive integer.");
    end

    settings.truncateFraction = 0.9;
    if isfield(S,"sawFFTTruncateFraction")
        settings.truncateFraction = double(S.sawFFTTruncateFraction);
    end
    if ~isscalar(settings.truncateFraction) || ...
            ~isfinite(settings.truncateFraction) || ...
            settings.truncateFraction <= 0 || ...
            settings.truncateFraction > 1
        error("TGS:SAWSettings", ...
            "S.sawFFTTruncateFraction must be in the interval (0, 1].");
    end

    settings.spectrumSmoothBins = 2*settings.zeroPadFactor+1;
    if isfield(S,"sawSpectrumSmoothBins")
        settings.spectrumSmoothBins = double(S.sawSpectrumSmoothBins);
    end
    if ~isscalar(settings.spectrumSmoothBins) || ...
            ~isfinite(settings.spectrumSmoothBins) || ...
            settings.spectrumSmoothBins < 1 || ...
            settings.spectrumSmoothBins ~= fix(settings.spectrumSmoothBins)
        error("TGS:SAWSettings", ...
            "S.sawSpectrumSmoothBins must be a positive integer.");
    end

    settings.minPeakRatio = 3.25;
    if isfield(S,"sawMinPeakRatio")
        settings.minPeakRatio = double(S.sawMinPeakRatio);
    end
    if ~isscalar(settings.minPeakRatio) || ...
            ~isfinite(settings.minPeakRatio) || settings.minPeakRatio <= 1
        error("TGS:SAWSettings", ...
            "S.sawMinPeakRatio must be a finite scalar greater than one.");
    end

    settings.singleTraceMinPeakRatio = 5;
    if isfield(S,"sawSingleTraceMinPeakRatio")
        settings.singleTraceMinPeakRatio = ...
            double(S.sawSingleTraceMinPeakRatio);
    end
    if ~isscalar(settings.singleTraceMinPeakRatio) || ...
            ~isfinite(settings.singleTraceMinPeakRatio) || ...
            settings.singleTraceMinPeakRatio <= 1
        error("TGS:SAWSettings", ...
            "S.sawSingleTraceMinPeakRatio must be greater than one.");
    end

    settings.bandFraction = 0.10;
    if isfield(S,"sawBandFraction")
        settings.bandFraction = double(S.sawBandFraction);
    end
    if ~isscalar(settings.bandFraction) || ...
            ~isfinite(settings.bandFraction) || ...
            settings.bandFraction <= 0 || settings.bandFraction >= 1
        error("TGS:SAWSettings", ...
            "S.sawBandFraction must be a finite scalar between zero and one.");
    end
    settings.minimumBandBins = 4;
    settings.minimumHalfBandwidthHz = ...
        settings.minimumBandBins/min(recordDuration);

    settings.halfBandwidthHz = NaN;
    if isfield(S,"sawHalfBandwidthHz") && ~isempty(S.sawHalfBandwidthHz)
        settings.halfBandwidthHz = double(S.sawHalfBandwidthHz);
        if ~isscalar(settings.halfBandwidthHz) || ...
                ~isfinite(settings.halfBandwidthHz) || ...
                settings.halfBandwidthHz <= 0
            error("TGS:SAWSettings", ...
                "S.sawHalfBandwidthHz must be a positive finite scalar.");
        end
    end

    settings.maskTransitionFraction = 0.25;
    if isfield(S,"sawMaskTransitionFraction")
        settings.maskTransitionFraction = ...
            double(S.sawMaskTransitionFraction);
    end
    if ~isscalar(settings.maskTransitionFraction) || ...
            ~isfinite(settings.maskTransitionFraction) || ...
            settings.maskTransitionFraction < 0 || ...
            settings.maskTransitionFraction >= 1
        error("TGS:SAWSettings", ...
            "S.sawMaskTransitionFraction must be in the interval [0, 1).");
    end

    settings.nativeResolutionHz = nativeResolution;
    settings.positiveLimitHz = positiveLimit;
    settings.configuredFrequencyHz = configuredFrequencyHz;
end

function [frequencyHz,peakRatio,requiredPeakRatio,detected] = ...
        detectSAWFrequency( ...
        sourceTime,sourceResidual,sampleInterval,settings)
% Follow the reference workflow: truncate, differentiate, FFT, and smooth.
% Repeat runs are normalized separately and then combined at one frequency.

    traceCount = numel(sourceTime);
    requiredPeakRatio = settings.minPeakRatio;
    if traceCount == 1
        requiredPeakRatio = max(requiredPeakRatio, ...
            settings.singleTraceMinPeakRatio);
    end
    paddedResolution = zeros(traceCount,1);
    frequencyCell = cell(traceCount,1);
    amplitudeCell = cell(traceCount,1);

    for traceIndex = 1:traceCount
        postPump = find(sourceTime{traceIndex} >= 0);
        analysisCount = min(numel(postPump),max(16, ...
            floor(settings.truncateFraction*numel(postPump))));
        if analysisCount < 16
            error("TGS:SAWAnalysisWindow", ...
                "FFT peak detection requires at least 16 post-pump samples.");
        end
        use = postPump(1:analysisCount);
        y = endpointDetrend(sourceResidual{traceIndex}(use));

        % Differentiation suppresses the slowly varying thermal background.
        derivative = diff(y)/sampleInterval(traceIndex);
        derivative = derivative-mean(derivative);
        sampleCount = numel(derivative);
        nfft = 2^nextpow2(settings.zeroPadFactor*sampleCount);
        powerSpectrum = abs(fft(derivative,nfft)).^2;
        positiveCount = floor((nfft-1)/2);
        frequency = (1:positiveCount).'/(nfft*sampleInterval(traceIndex));
        powerSpectrum = powerSpectrum(2:positiveCount+1);

        smoothBins = min(settings.spectrumSmoothBins,positiveCount);
        if mod(smoothBins,2) == 0 && smoothBins > 1
            smoothBins = smoothBins-1;
        end
        if smoothBins > 1
            powerSpectrum = movmean(powerSpectrum,smoothBins, ...
                "Endpoints","shrink");
        end
        amplitude = sqrt(max(powerSpectrum,0));
        searchUse = frequency >= settings.searchBandHz(1) & ...
            frequency <= settings.searchBandHz(2);
        if nnz(searchUse) < 5
            error("TGS:SAWSearchBand", ...
                "The SAW search band contains fewer than five FFT samples.");
        end
        scale = median(amplitude(searchUse));
        if ~isfinite(scale) || scale <= 0
            scale = max(amplitude(searchUse));
        end
        if ~isfinite(scale) || scale <= 0
            scale = 1;
        end
        frequencyCell{traceIndex} = frequency;
        amplitudeCell{traceIndex} = amplitude/scale;
        paddedResolution(traceIndex) = ...
            1/(nfft*sampleInterval(traceIndex));
    end

    gridStep = min(paddedResolution);
    firstCommonFrequency = max(cellfun(@(frequency) ...
        frequency(1),frequencyCell));
    lastCommonFrequency = min(cellfun(@(frequency) ...
        frequency(end),frequencyCell));
    gridStart = max(settings.searchBandHz(1),firstCommonFrequency);
    gridEnd = min(settings.searchBandHz(2),lastCommonFrequency);
    frequencyGrid = (gridStart:gridStep:gridEnd).';
    if numel(frequencyGrid) < 5
        error("TGS:SAWSearchBand", ...
            "The shared SAW search band contains fewer than five FFT samples.");
    end

    normalizedAmplitude = zeros(numel(frequencyGrid),traceCount);
    for traceIndex = 1:traceCount
        normalizedAmplitude(:,traceIndex) = interp1( ...
            frequencyCell{traceIndex},amplitudeCell{traceIndex}, ...
            frequencyGrid,"linear");
    end
    if any(~isfinite(normalizedAmplitude),"all")
        error("TGS:SAWSpectrum", ...
            "The normalized FFT spectra could not be placed on a common grid.");
    end
    aggregateAmplitude = median(normalizedAmplitude,2);

    localPeak = find(aggregateAmplitude(2:end-1) > ...
        aggregateAmplitude(1:end-2) & ...
        aggregateAmplitude(2:end-1) >= ...
        aggregateAmplitude(3:end))+1;
    if isempty(localPeak)
        frequencyHz = NaN;
        peakRatio = 0;
        detected = false;
        return
    end

    edgeGuardHz = max(settings.minimumHalfBandwidthHz, ...
        settings.bandFraction*frequencyGrid(localPeak));
    awayFromEdge = frequencyGrid(localPeak)-settings.searchBandHz(1) > ...
        edgeGuardHz & settings.searchBandHz(2)-frequencyGrid(localPeak) > ...
        edgeGuardHz;
    localPeak = localPeak(awayFromEdge);
    if isempty(localPeak)
        frequencyHz = NaN;
        peakRatio = 0;
        detected = false;
        return
    end

    [peakAmplitude,selected] = max(aggregateAmplitude(localPeak));
    peakIndex = localPeak(selected);
    spectralBaseline = median(aggregateAmplitude);
    if ~isfinite(spectralBaseline) || spectralBaseline <= 0
        frequencyHz = NaN;
        peakRatio = 0;
        detected = false;
        return
    end

    % The FFT-bin maximum is used directly; no acoustic peak is fitted.
    frequencyHz = frequencyGrid(peakIndex);
    peakRatio = peakAmplitude/spectralBaseline;
    detected = peakRatio >= requiredPeakRatio;
end

function [sawTerm,powerBefore,powerAfter] = extractSAWBand( ...
        t,residual,dt,frequencyHz,halfBandwidthHz,settings)
% Select conjugate frequency bands and transform them back to the time domain.

    postPump = find(t >= 0);
    if numel(postPump) < 2
        error("TGS:SAWAnalysisWindow", ...
            "IFFT reconstruction requires at least two post-pump samples.");
    end

    postResidual = endpointDetrend(residual(postPump));
    sampleCount = numel(postResidual);
    nfft = 2^nextpow2(3*sampleCount);
    firstSample = sampleCount;
    paddedResidual = zeros(nfft,1);
    paddedResidual(firstSample:firstSample+sampleCount-1) = postResidual;
    sampleRate = 1/dt;
    frequency = (0:nfft-1).'*(sampleRate/nfft);
    frequency(frequency > sampleRate/2) = ...
        frequency(frequency > sampleRate/2)-sampleRate;
    distance = abs(abs(frequency)-frequencyHz);

    transitionWidth = ...
        settings.maskTransitionFraction*halfBandwidthHz;
    coreWidth = halfBandwidthHz-transitionWidth;
    mask = zeros(nfft,1);
    mask(distance <= coreWidth) = 1;
    if transitionWidth > 0
        transition = distance > coreWidth & ...
            distance < halfBandwidthHz;
        mask(transition) = 0.5*(1+cos(pi* ...
            (distance(transition)-coreWidth)/transitionWidth));
    end

    if ~any(mask > 0)
        error("TGS:SAWBandwidth", ...
            "The SAW mask does not contain a native FFT bin.");
    end

    spectrum = fft(paddedResidual);
    paddedTerm = ifft(spectrum.*mask,"symmetric");
    sawTerm = zeros(size(residual));
    sawTerm(postPump) = real(paddedTerm( ...
        firstSample:firstSample+sampleCount-1));

    % Zero padding makes the time-domain convolution non-circular over the
    % measured interval. Pre-pump samples are never altered.
    if any(~isfinite(sawTerm))
        error("TGS:SAWTerm", ...
            "The inverse-FFT SAW term contains nonfinite values.");
    end

    powerBefore = bandPower(t,residual,dt,frequencyHz,halfBandwidthHz);
    powerAfter = bandPower(t,residual-sawTerm,dt, ...
        frequencyHz,halfBandwidthHz);
end

function power = bandPower(t,y,dt,frequencyHz,halfBandwidthHz)
% Measure post-pump FFT power in the reported removal band.

    use = t >= 0;
    y = endpointDetrend(y(use));
    sampleCount = numel(y);
    if sampleCount < 2
        power = NaN;
        return
    end
    sampleRate = 1/dt;
    frequency = (0:sampleCount-1).'*(sampleRate/sampleCount);
    positiveBand = frequency >= frequencyHz-halfBandwidthHz & ...
        frequency <= frequencyHz+halfBandwidthHz;
    spectrum = fft(y);
    power = sum(abs(spectrum(positiveBand)).^2)/sampleCount^2;
end

function y = endpointDetrend(y)
% Remove the endpoint-to-endpoint line to limit spectral leakage.

    y = y(:);
    if numel(y) < 2
        return
    end
    position = linspace(0,1,numel(y)).';
    endpointLine = y(1)+(y(end)-y(1))*position;
    y = y-endpointLine;
end


