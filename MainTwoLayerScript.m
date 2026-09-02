function results = MainTwoLayerScript()
% Fit alpha_f, alpha_s, and R for one spot and grating
% Edit the user-input section, then call results = MainTwoLayerScript

%% User inputs: files and naming convention
    S.dataDir = "C:\Users\Maken\OneDrive - Massachusetts Institute of Technology\FFUSars\MultilayerV2-model-work\MultilayerV2-model-work\tungsten_depth_study";
    S.calFile = "C:\Users\Maken\OneDrive - Massachusetts Institute of Technology\FFUSars\MultilayerV2-model-work\MultilayerV2-model-work\tungsten_depth_study\Tungsten_Calibration_06.40um_postprocessing.txt";
    S.filePattern = "*.txt";
    S.fileRegex = "(?<nominal>\d+\.\d+)(?:um)?[_-]" + ...
        "(?<location>baseline|spot\d+(?:-baseline)?)-" + ...
        "(?<polarity>POS|NEG)-(?<run>\d+)\.txt$";

%% User inputs: measurement-file layout
    S.nHeader = 16;       % Header lines skipped in each measurement file
    S.tCol = 1;           % Time column
    S.yCol = 2;           % Detector-signal column
    S.tScale = 1;         % Conversion from file time units to seconds

%% User inputs: pump arrival and fit window
    S.tSearch = [0,100e-9];  % Pump-arrival search interval [s]
    S.nShift = 0;            % Samples retained before detected arrival
    S.tFit = [0,inf];         % Preprocessing interval relative to arrival [s]
    S.nFitPoints = 800;       % Maximum nonlinear-fit samples per run
    S.smoothTime = 4e-9;      % Thermal-prefit moving-average width [s]
    S.fitStartTime = 0;       % Additional final-fit delay [s]

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

%% User inputs: mandatory [alpha_f, alpha_s, R] fit
    % Every fit always optimizes all three values in this order.
    S.p0 = [6.7e-5,2.7e-5,1e-10];
    S.pLower = [1.5e-7,1.5e-7,1e-14];
    S.pUpper = [1e-1,1e-1,1e-6];

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
    S.boundTolerance = 1e-6; % Bound test in log10 parameter coordinates
    S.dx = 2e-3;             % Finite-difference step in log10 coordinates
    S.display = "off";       % Solver output; run progress still prints

%% User inputs: plotting
    S.makePlots = true;
    S.figureVisible = "on";

%% Load the one selected spot and grating
    [traces,spot] = prepareTGSData(S);
    if isempty(traces)
        error("TGS:Data", ...
            "No complete calibrated POS/NEG/baseline sets were found.");
    end

    if any([traces.spot] ~= spot)
        error("TGS:SpotSelection", ...
            "Prepared runs do not match the calibration-selected spot.");
    end

%% Fit each run and combine repeated estimates
    fit = fitTwoLayerTGSRuns(traces,S);

%% Package the spot-level result
    results.spot = spot;
    results.alpha_f = fit.p(1);
    results.alpha_s = fit.p(2);
    results.R = fit.p(3);
    results.x = fit.x;
    results.sensitivity = fit.sensitivity;
    results.parameterError = fit.parameterError;
    results.withinFitError = fit.withinFitError;
    results.betweenRunStd = fit.betweenRunStd;
    results.parameterRunCount = fit.parameterRunCount;
    results.parameterIdentifiable = fit.parameterIdentifiable;
    results.identifiable = fit.identifiable;
    results.traces = fit.traces;
    results.runFits = fit.fits;
    results.runSummary = fit.runSummary;
    results.runCount = fit.runCount;
    results.r = fit.r;
    results.resnorm = fit.resnorm;
    results.objectiveResnorm = fit.objectiveResnorm;
    results.exitflag = fit.exitflag;
    results.output = fit.output;
    results.atLowerBound = fit.atLowerBound;
    results.atUpperBound = fit.atUpperBound;
    results.jacobianRank = fit.jacobianRank;
    results.jacobianCondition = fit.jacobianCondition;
    results.degreesOfFreedom = fit.degreesOfFreedom;
    results.physicalJacobianRank = fit.physicalJacobianRank;
    results.physicalJacobianCondition = fit.physicalJacobianCondition;
    results.physicalDegreesOfFreedom = fit.physicalDegreesOfFreedom;

%% Report parameter and run diagnostics
    parameterName = ["alpha_f";"alpha_s";"R"];
    unit = ["m^2/s";"m^2/s";"m^2 K/W"];
    results.summary = table( ...
        repmat(spot,3,1),parameterName,fit.p(:),fit.sensitivity(:), ...
        fit.parameterError(:),fit.withinFitError(:), ...
        fit.betweenRunStd(:),fit.parameterRunCount(:), ...
        repmat(fit.runCount,3,1),fit.parameterIdentifiable(:), ...
        fit.atLowerBound(:),fit.atUpperBound(:),unit, ...
        'VariableNames',{'Spot','Symbol','Value','LocalSensitivity', ...
        'CombinedError','WithinFitError','BetweenRunStd','RunsUsed', ...
        'RunsTotal','Identifiable','AnyAtLowerBound', ...
        'AnyAtUpperBound','Unit'});

    fprintf("spot %d, %d run(s), 1 grating\n",spot,fit.runCount);
    disp(results.summary)
    disp(results.runSummary)
    fprintf("successful optimizer exits = %d/%d\n", ...
        nnz(fit.exitflag > 0),fit.runCount);
    fprintf("unweighted ||r||_2 = %.6g\n",norm(fit.r));
    fprintf("full-rank physical Jacobians = %d/%d\n", ...
        nnz(fit.physicalJacobianRank == 3),fit.runCount);

%% Plot run-level fits and physical sensitivities
    if S.makePlots
        results.figure = plotTwoLayerTGSResults(results,S);
    else
        results.figure = gobjects(0);
    end
end
