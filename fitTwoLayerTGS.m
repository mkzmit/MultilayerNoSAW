function fit = fitTwoLayerTGS(trace,S,initializationTrace)
% Fit alpha_f, alpha_s, and R for one run.
% Inputs:
%   trace - Scalar prepared-trace structure from prepareTGSData
%   S - Model, bounds, FFT-initialization, and optimizer settings
%   initializationTrace - Optional FFT-smoothed copy of trace. Repeated-run
%                         callers use this to share one detected SAW peak.
% Output:
%   fit - Thermal-model parameters and selected multistart trial index.

%% Validate 
    if ~isscalar(trace)
        error("TGS:OneRunFit", "fitTwoLayerTGS fits one run at a time; use fitTwoLayerTGSRuns for repeats.");
    end
    [x0,lowerX,upperX] = physicalSearchSpace(S);
    options = fittingOptions(S);

%% Build an FFT-smoothed trace only for physical-parameter initialization
    if nargin < 3 || isempty(initializationTrace)
        initializationTrace = fftThermalInitialization(trace,S);
    else
        if ~isscalar(initializationTrace) || ...
                ~isfield(initializationTrace,"fftInitialized") || ...
                ~isequal(initializationTrace.fftInitialized,true) || ...
                ~isfield(initializationTrace,"fftInitialization")
            error("TGS:InitializationTrace", ...
                "The supplied initialization trace must come from fftThermalInitialization.");
        end
    end

    if numel(initializationTrace.t) ~= numel(trace.t) || ...
            any(initializationTrace.t(:) ~= trace.t(:))
        error("TGS:InitializationTrace", ...
            "The FFT initialization and final-fit traces must use the same fit samples.");
    end

    initializationWeight = traceWeight(initializationTrace);
    initializationData = initializationWeight*initializationTrace.y(:);
    dummyData = zeros(size(initializationData));
    initializationModel = @(x,xdata) weightedThermalModel( x,xdata,initializationTrace,S,initializationWeight);
    try
        xThermal = lsqcurvefit(initializationModel,x0,dummyData, ...
            initializationData,lowerX,upperX,options);
    catch exception
        if strcmp(exception.identifier,"MATLAB:OperationTerminated")
            rethrow(exception)
        end
        warning("TGS:Initialization", ...
            "Thermal prefit failed (%s); continuing with multistart points.", ...
            exception.message);
        xThermal = x0;
    end
    xThermal = xThermal(:).';

%% Fit the original unsmoothed response without an added SAW term
    fitTrace = finalFitTrace(trace,S);
    finalStarts = multistartPoints([xThermal;x0],lowerX,upperX,S);
    finalWeight = traceWeight(fitTrace);
    measuredWeighted = finalWeight*fitTrace.y(:);
    dummyData = zeros(size(measuredWeighted));
    thermalModel = @(x,xdata) weightedThermalModel( x,xdata,fitTrace,S,finalWeight);

    startResults = repmat(struct("x",nan(1,3), ...
        "resnorm",Inf,"exitflag",NaN,"errorMessage",""), ...
        size(finalStarts,1),1);
    bestResnorm = Inf;
    selectedStart = NaN;
    bestConverged = false;

    for startIndex = 1:size(finalStarts,1)
        try
            [candidateX,candidateResnorm,~,candidateExitflag] = ...
                lsqcurvefit(thermalModel,finalStarts(startIndex,:), ...
                dummyData,measuredWeighted,lowerX,upperX,options);
        catch exception
            if strcmp(exception.identifier,"MATLAB:OperationTerminated")
                rethrow(exception)
            end
            startResults(startIndex).errorMessage = string(exception.message);
            continue
        end

        startResults(startIndex).x = candidateX(:).';
        startResults(startIndex).resnorm = candidateResnorm;
        startResults(startIndex).exitflag = candidateExitflag;

        candidateIsFinite = isfinite(candidateResnorm) && ...
            all(isfinite(candidateX));
        candidateConverged = candidateIsFinite && candidateExitflag > 0;
        candidateIsBetter = candidateIsFinite && ...
            (~isfinite(selectedStart) || ...
            (candidateConverged && ~bestConverged) || ...
            (candidateConverged == bestConverged && ...
            candidateResnorm < bestResnorm));

        if candidateIsBetter
            bestX = candidateX(:).';
            bestResnorm = candidateResnorm;
            selectedStart = startIndex;
            bestConverged = candidateConverged;
        end
    end

    if ~isfinite(selectedStart)
        errorMessages = string({startResults.errorMessage});
        firstError = find(strlength(errorMessages) > 0,1);
        if isempty(firstError)
            detail = "No trial returned finite parameters and an objective.";
        else
            detail = "First trial error: " + errorMessages(firstError);
        end
        error("TGS:Optimization", ...
            "All %d multistart trials failed. %s",numel(startResults),detail);
    end

    successfulStart = arrayfun(@(trial) trial.exitflag > 0 && ...
        isfinite(trial.resnorm) && all(isfinite(trial.x)),startResults);
    if ~any(successfulStart)
        warning("TGS:Optimization", ...
            "No multistart trial reported convergence; using the lowest finite objective.");
    end

%% Package the per-run result
    fit.p = 10.^bestX;
    fit.selectedStart = selectedStart;
end

function [x0,lowerX,upperX] = physicalSearchSpace(S)
% Validate and log-transform the three physical inputs.

    requiredFields = ["p0","pLower","pUpper"];
    for fieldName = requiredFields
        if ~isfield(S,fieldName)
            error("TGS:Parameters","Missing required setting S.%s.",fieldName);
        end
    end

    p0 = double(S.p0(:).');
    lower = double(S.pLower(:).');
    upper = double(S.pUpper(:).');
    if numel(p0) ~= 3 || numel(lower) ~= 3 || numel(upper) ~= 3
        error("TGS:Parameters", ...
            "S.p0, S.pLower, and S.pUpper must each contain [alpha_f, alpha_s, R].");
    end

    if any(~isfinite([p0,lower,upper])) || any([p0,lower,upper] <= 0)
        error("TGS:Parameters", ...
            "All physical initial values and bounds must be positive and finite.");
    end

    if any(lower >= upper) || any(p0 <= lower) || any(p0 >= upper)
        error("TGS:Parameters", ...
            "Each initial physical value must lie strictly inside its lower and upper bounds.");
    end

    x0 = log10(p0);
    lowerX = log10(lower);
    upperX = log10(upper);
end

function starts = multistartPoints(primaryStarts,lowerX,upperX,S)
% Build reproducible, space-filling starts in log10 parameter coordinates.

    startCount = 24;
    if isfield(S,"multistartCount")
        startCount = double(S.multistartCount);
    end
    if ~isscalar(startCount) || ~isfinite(startCount) || ...
            startCount < 4 || startCount ~= fix(startCount)
        error("TGS:Multistart", ...
            "S.multistartCount must be an integer of at least 4.");
    end

    seed = 1;
    if isfield(S,"multistartSeed")
        seed = double(S.multistartSeed);
    end
    if ~isscalar(seed) || ~isfinite(seed) || seed < 0 || ...
            seed > double(intmax("uint32")) || seed ~= fix(seed)
        error("TGS:Multistart", ...
            "S.multistartSeed must be an integer from 0 through 2^32-1.");
    end

    validPrimary = all(isfinite(primaryStarts),2) & ...
        all(primaryStarts >= lowerX & primaryStarts <= upperX,2);
    midpoint = (lowerX+upperX)/2;
    anchors = unique([primaryStarts(validPrimary,:);midpoint], ...
        "rows","stable");
    anchors = anchors(1:min(startCount,size(anchors,1)),:);
    spaceCount = startCount-size(anchors,1);

    if spaceCount == 0
        starts = anchors;
        return
    end

    % Stratify every coordinate independently, then permute its strata.
    stream = RandStream("mt19937ar","Seed",seed);
    unitPoints = zeros(spaceCount,numel(lowerX));
    for parameterIndex = 1:numel(lowerX)
        strata = randperm(stream,spaceCount).';
        unitPoints(:,parameterIndex) = ...
            (strata-1+rand(stream,spaceCount,1))/spaceCount;
    end

    spaceStarts = lowerX+unitPoints.*(upperX-lowerX);
    starts = [anchors;spaceStarts];
end

function fitTrace = finalFitTrace(trace,S)
% Apply an optional fixed post-pump delay to the final fit.

    startTime = 0;
    if isfield(S,"fitStartTime")
        startTime = S.fitStartTime;
    end
    if ~isscalar(startTime) || ~isfinite(startTime)
        error("TGS:FitWindow","S.fitStartTime must be a finite scalar.");
    end

    t = trace.t(:);
    y = trace.y(:);
    use = t >= startTime;
    if nnz(use) < 10
        error("TGS:FitWindow", "The final fit window contains fewer than 10 points.");
    end

    fitTrace = trace;
    fitTrace.t = t(use);
    fitTrace.y = y(use);
end

function options = fittingOptions(S)
% Build bounded least-squares options for both thermal stages

    options = optimoptions("lsqcurvefit", ...
        "Algorithm","trust-region-reflective", ...
        "FiniteDifferenceType","forward", ...
        "FiniteDifferenceStepSize",S.dx, ...
        "MaxIterations",S.maxIter, ...
        "MaxFunctionEvaluations",S.maxEval, ...
        "FunctionTolerance",S.ftol, ...
        "StepTolerance",S.xtol, ...
        "OptimalityTolerance",S.gtol, ...
        "Display",S.display);
end

function weightedFit = weightedThermalModel(x,~,trace,S,weight)
% Profile temperature, displacement, and offset terms

    p = physicalParameters(x);
    [temperature,displacement] = TwoLayerModel(trace.t,trace.Lambda,p,S);
    design = [temperature,displacement,ones(size(trace.t(:)))];
    yfit = linearFit(design,trace.y);

    weightedFit = weight*yfit;
end

function p = physicalParameters(x)
% Convert three log10 coordinates to physical values.

    x = x(:).';
    if numel(x) ~= 3 || any(~isfinite(x))
        error("TGS:Parameters", ...
            "The optimizer must supply log10([alpha_f, alpha_s, R]).");
    end
    p = 10.^x;
end

function yfit = linearFit(design,y)
% Solve scaled linear nuisance amplitudes at fixed physical values.

    y = y(:);
    if size(design,1) ~= numel(y)
        error("TGS:ModelSize", ...
            "The design matrix and measured run have different lengths.");
    end
    if any(~isfinite(design),"all") || any(~isfinite(y))
        error("TGS:Model", ...
            "The model response and measured data must be finite.");
    end

    scale = vecnorm(design,2,1);
    if any(~isfinite(scale)) || any(scale == 0)
        error("TGS:Model", ...
            "At least one model component has zero or nonfinite magnitude.");
    end

    scaledDesign = design./scale;
    scaledCoefficients = scaledDesign\y;
    coefficients = scaledCoefficients./scale.';
    yfit = design*coefficients;
end

function weight = traceWeight(trace)
% Scale one residual by measured pre-pump noise when available

    y = trace.y(:);
    if any(~isfinite(y))
        error("TGS:DataScale", "The measured run contains nonfinite values.");
    end

    scale = NaN;
    
    if isfield(trace,"noiseStd") && isscalar(trace.noiseStd) && ...
            isfinite(trace.noiseStd) && trace.noiseStd > 0
        scale = trace.noiseStd;
    end
    
    if ~isfinite(scale) || scale <= 0
        scale = norm(y-mean(y));
    end
    
    if scale == 0
        error("TGS:DataScale", "The measured run has no signal variation.");
    end
    
    weight = 1/scale;
end
