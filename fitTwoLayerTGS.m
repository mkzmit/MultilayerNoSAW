function fit = fitTwoLayerTGS(trace,S)
% Fit alpha_f, alpha_s, and R for one  run
% Inputs:
%   trace - Scalar prepared-trace structure from prepareTGSData
%   S - Model, bounds, smoothing, and optimizer settings
% Output:
%   fit - Thermal-model parameters, reconstructed signal, and diagnostics

%% Validate 
    if ~isscalar(trace)
        error("TGS:OneRunFit", "fitTwoLayerTGS fits one run at a time; use fitTwoLayerTGSRuns for repeats.");
    end
    [x0,lowerX,upperX] = physicalSearchSpace(S);
    options = fittingOptions(S);

%% Fit a smoothed trace to initialize the physical parameters
    initializationTrace = smoothThermalTrace(trace,S);
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

%% Fit the unsmoothed thermal and displacement response
    fitTrace = finalFitTrace(trace,S);
    finalStarts = multistartPoints([xThermal;x0],lowerX,upperX,S);
    finalWeight = traceWeight(fitTrace);
    measuredWeighted = finalWeight*fitTrace.y(:);
    dummyData = zeros(size(measuredWeighted));
    thermalModel = @(x,xdata) weightedThermalModel( x,xdata,fitTrace,S,finalWeight);

    startResults = repmat(struct("x0",[],"x",nan(1,3), ...
        "resnorm",Inf,"exitflag",NaN,"output",[], ...
        "errorIdentifier","","errorMessage",""),size(finalStarts,1),1);
    bestResnorm = Inf;
    selectedStart = NaN;
    bestConverged = false;

    for startIndex = 1:size(finalStarts,1)
        startResults(startIndex).x0 = finalStarts(startIndex,:);
        try
            [candidateX,candidateResnorm,candidateResidual, ...
                candidateExitflag,candidateOutput,~,candidateJacobian] = ...
                lsqcurvefit(thermalModel,finalStarts(startIndex,:), ...
                dummyData,measuredWeighted,lowerX,upperX,options);
        catch exception
            if strcmp(exception.identifier,"MATLAB:OperationTerminated")
                rethrow(exception)
            end
            startResults(startIndex).errorIdentifier = string(exception.identifier);
            startResults(startIndex).errorMessage = string(exception.message);
            continue
        end

        startResults(startIndex).x = candidateX(:).';
        startResults(startIndex).resnorm = candidateResnorm;
        startResults(startIndex).exitflag = candidateExitflag;
        startResults(startIndex).output = candidateOutput;

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
            objectiveResidual = candidateResidual;
            exitflag = candidateExitflag;
            output = candidateOutput;
            jacobian = candidateJacobian;
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
    successfulStartCount = nnz(successfulStart);
    if successfulStartCount == 0
        warning("TGS:Optimization", ...
            "No multistart trial reported convergence; using the lowest finite objective.");
    end

%% Reconstruct the unweighted thermal-model signal
    [yfit,rebuilt] = weightedThermalModel(bestX,[],fitTrace,S,1);
    measured = fitTrace.y(:);
    physicalValues = 10.^bestX;

%% Calculate uncertainty and identifiability diagnostics
    boundTolerance = S.dx;
    if isfield(S,"boundTolerance") && isfinite(S.boundTolerance) && S.boundTolerance >= 0
        boundTolerance = S.boundTolerance;
    end

    atLowerBound = abs(bestX-lowerX) <= boundTolerance;
    atUpperBound = abs(bestX-upperX) <= boundTolerance;
    physicalAtBound = atLowerBound | atUpperBound;
    statistics = nonlinearStatistics( ...
        jacobian,bestResnorm,rebuilt.linearRank);
    xError = statistics.xError;
    xError(physicalAtBound) = NaN;
    statisticsValid = statistics.fullRank && ...
        statistics.degreesOfFreedom > 0;
    parameterIdentifiable = statisticsValid & ~physicalAtBound;

%% Package the per-run result
    fit.x = bestX;
    fit.p = physicalValues;
    fit.optimizerX = bestX;
    fit.resnorm = bestResnorm;
    fit.r = measured-yfit;
    fit.objectiveResidual = objectiveResidual;
    fit.yfit = yfit;
    fit.traces = rebuilt;
    fit.exitflag = exitflag;
    fit.output = output;
    fit.startResults = startResults;
    fit.selectedStart = selectedStart;
    fit.successfulStartCount = successfulStartCount;
    fit.J = jacobian;
    fit.atLowerBound = atLowerBound;
    fit.atUpperBound = atUpperBound;
    fit.parameterError = log(10)*physicalValues.*xError.';
    fit.sensitivity = vecnorm(jacobian,2,1);
    fit.parameterIdentifiable = parameterIdentifiable;
    fit.identifiable = all(parameterIdentifiable);
    fit.jacobianColumnNorm = statistics.columnNorm;
    fit.jacobianSingularValues = statistics.singularValues;
    fit.jacobianRank = statistics.rank;
    fit.jacobianCondition = statistics.condition;
    fit.degreesOfFreedom = statistics.degreesOfFreedom;
    fit.covariance = statistics.covariance;
    fit.correlation = statistics.correlation;
    fit.physicalJacobianSingularValues = statistics.singularValues;
    fit.physicalJacobianRank = statistics.rank;
    fit.physicalJacobianCondition = statistics.condition;
    fit.physicalDegreesOfFreedom = statistics.degreesOfFreedom;
    fit.physicalCovariance = statistics.covariance;
    fit.physicalCorrelation = statistics.correlation;
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

function thermalTrace = smoothThermalTrace(trace,S)
% Smooth the full-resolution signal for initialization

    if isfield(trace,"tFull") && isfield(trace,"yFull") &&  ~isempty(trace.tFull) && ~isempty(trace.yFull)
        sourceTime = trace.tFull(:);
        sourceSignal = trace.yFull(:);
    else
        sourceTime = trace.t(:);
        sourceSignal = trace.y(:);
    end

    timeSteps = diff(sourceTime);
    
    if numel(sourceTime) < 3 || any(~isfinite(sourceTime)) || any(~isfinite(sourceSignal)) || any(timeSteps <= 0)
        error("TGS:Trace", "The run must contain finite, increasing time samples.");
    end
    
    if ~isscalar(S.smoothTime) || ~isfinite(S.smoothTime) || S.smoothTime < 0
        error("TGS:Smoothing", "S.smoothTime must be a nonnegative finite duration.");
    end

    if S.smoothTime == 0
        smoothedSignal = sourceSignal;
    else
        smoothedSignal = movmean(sourceSignal,S.smoothTime,"SamplePoints",sourceTime,"Endpoints","shrink");
    end

    thermalTrace = trace;
    thermalTrace.t = trace.t(:);
    thermalTrace.y = interp1(sourceTime,smoothedSignal,thermalTrace.t,"linear");
end

function [weightedFit,state] = weightedThermalModel( x,~,trace,S,weight)
% Profile temperature, displacement, and offset terms

    p = physicalParameters(x);
    [temperature,displacement] = TwoLayerModel(trace.t,trace.Lambda,p,S);
    design = [temperature,displacement,ones(size(trace.t(:)))];
    [yfit,coefficients] = linearFit(design,trace.y);

    weightedFit = weight*yfit;
    state = trace;
    state.Theta = temperature;
    state.Uz = displacement;
    state.c = coefficients;
    state.yfit = yfit;
    state.r = trace.y(:)-yfit;
    state.linearRank = rank(design);
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

function [yfit,coefficients] = linearFit(design,y)
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

function statistics = nonlinearStatistics(J,resnorm,linearParameterCount)
% Calculate covariance and scaled-Jacobian rank metrics

    J = full(J);
    parameterCount = size(J,2);
    columnNorm = vecnorm(J,2,1);
    validColumn = isfinite(columnNorm) & columnNorm > 0;
    scaledJ = zeros(size(J));
    scaledJ(:,validColumn) = J(:,validColumn)./columnNorm(validColumn);
    [~,singularMatrix,rightVectors] = svd(scaledJ,"econ");
    singularValues = diag(singularMatrix);

    if isempty(singularValues) || ~isfinite(singularValues(1))
        nonlinearRank = 0;
        conditionNumber = Inf;
    else
        tolerance = max(size(scaledJ))*eps(singularValues(1));
        nonlinearRank = nnz(singularValues > tolerance);
        if nonlinearRank == parameterCount && singularValues(end) > 0
            conditionNumber = singularValues(1)/singularValues(end);
        else
            conditionNumber = Inf;
        end
    end

    degreesOfFreedom = size(J,1)-linearParameterCount-nonlinearRank;
    fullRank = all(validColumn) && nonlinearRank == parameterCount;
    covariance = nan(parameterCount);
    correlation = nan(parameterCount);
    xError = nan(parameterCount,1);

    if fullRank && degreesOfFreedom > 0
        inverseScale = diag(1./columnNorm);
        covariance = (resnorm/degreesOfFreedom)*(inverseScale* ...
            rightVectors*diag(1./singularValues.^2)*rightVectors.'* ...
            inverseScale);
        covariance = (covariance+covariance.')/2;
        xError = sqrt(max(diag(covariance),0));
        standardDeviation = sqrt(diag(covariance));
        correlation = covariance./(standardDeviation*standardDeviation.');
    end

    statistics.columnNorm = columnNorm;
    statistics.singularValues = singularValues;
    statistics.rank = nonlinearRank;
    statistics.condition = conditionNumber;
    statistics.degreesOfFreedom = degreesOfFreedom;
    statistics.fullRank = fullRank;
    statistics.covariance = covariance;
    statistics.correlation = correlation;
    statistics.xError = xError;
end
