function results = MainTwoLayerScript()
% Fit alpha_f, alpha_s, and R for each run from one spot and grating.
% Edit the user-input section, then call results = MainTwoLayerScript.

%% User inputs: files and naming convention
    S.dataDir = "C:\Users\Maken\OneDrive - Massachusetts Institute of Technology\FFUSars\TGS Data rep-20261006T182917Z-1-001\TGS Data rep\MKZ_coatings\WonNb\2026-08-11";
    S.calFile = "C:\Users\Maken\OneDrive - Massachusetts Institute of Technology\FFUSars\TGS Data rep-20261006T182917Z-1-001\TGS Data rep\MKZ_coatings\WonNb\2026-08-11\WonNb_spot1_postprocessing.txt";
    S.filePattern = "*.txt";
    S.fileRegex = "(?<nominal>\d+\.\d+)(?:um)?[_-]" + "(?<location>baseline|spot\d+(?:-baseline)?)-" + "(?<polarity>POS|NEG)-(?<run>\d+)\.txt$";

%% User inputs: measurement-file layout
    S.nHeader = 16;       % Header lines skipped in each measurement file
    S.tCol = 1;           % Time column
    S.yCol = 2;           % Detector-signal column
    S.tScale = 1;         % Conversion from file time units to seconds

%% User inputs: pump arrival and fit window
    S.tSearch = [0,100e-9];  % Pump-arrival search interval [s]
    S.nShift = 0;            % Samples retained before detected arrival
    S.tFit = [0,inf];        % Preprocessing interval relative to arrival [s]
    S.nFitPoints = 800;      % Maximum nonlinear-fit samples per run
    S.fitStartTime = 0;      % Additional final-fit delay [s]

%% User inputs: FFT smoothing used only for thermal initialization
    S.sawVelocityRange = [1.5e3,3e3]; % SAW search range [m/s]
    S.sawBandFraction = 0.10;          % IFFT half-band / detected frequency
    S.sawMinPeakRatio = 2.0;           % Repeated-run peak / median amplitude
    S.sawFFTTruncateFraction = 0.9;    % Post-pump fraction used for detection

%% User inputs: known film and substrate properties
    S.L = 1.5e-6;           % Film thickness [m]
    S.Cf = 2.56e6;          % Film volumetric heat capacity [J m^-3 K^-1]
    S.Cs = 2.27e6;          % Substrate volumetric heat capacity [J m^-3 K^-1]
    S.Ef = 4.11e11;         % Film Young's modulus [Pa]
    S.Es = 1.05e11;         % Substrate Young's modulus [Pa]
    S.nuf = 0.28;           % Film Poisson ratio
    S.nus = 0.40;           % Substrate Poisson ratio
    S.alphathf = 4.42e-6;   % Film thermal expansion coefficient [K^-1]
    S.alphaths = 7.07e-6;   % Substrate thermal expansion coefficient [K^-1]

%% User inputs: [alpha_f, alpha_s, R] guess and bounds
    S.p0 = [6.7e-5,2.7e-5,1e-8];
    S.pLower = [1.5e-7,1.5e-7,1e-14];
    S.pUpper = [1e-4,1e-4,1e-6];

%% User inputs: Fourier reconstruction
    S.Nw = 256;             % Even number of angular-frequency samples
    S.wmin = 2*pi*1e3;      % Minimum angular frequency [rad/s]
    S.wmax = 2*pi*1e12;     % Maximum angular frequency [rad/s]
    S.Q0 = 1;               % Unit absorbed surface impulse

%% User inputs: bounded nonlinear least squares
    S.maxIter = 1e10;
    S.maxEval = 4e10;
    S.ftol = 1e-8;
    S.xtol = 1e-8;
    S.gtol = 1e-8;
    S.dx = 2e-3;             % Finite-difference step in log10 coordinates
    S.multistartCount = 24;   % Total final-fit starts per run
    S.multistartSeed = 1;     % Reproducible log-space start-point sampling
    S.display = "off";        % Solver output

%% Load the one selected spot and grating
    [traces,spot] = prepareTGSData(S);
    if isempty(traces)
        error("TGS:Data", "No complete calibrated POS/NEG/baseline sets were found.");
    end

    if any([traces.spot] ~= spot)
        error("TGS:SpotSelection", "Prepared runs do not match the calibration-selected spot.");
    end

%% Fit each run independently
    results = fitTwoLayerTGSRuns(traces,S);

%% Report only the requested per-run values
    disp(results)
end
