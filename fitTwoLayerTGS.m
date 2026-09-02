function fit = fitTwoLayerTGS(trace,S)
%FITTWOLAYERTGS Fit alpha_f, alpha_s, and R for one measured run.
% Inputs:
%   trace - Scalar prepared-trace structure from prepareTGSData.
%   S     - Model, bounds, smoothing, and optimizer settings.
% Output:
%   fit   - Thermal-model parameters, reconstructed signal, and diagnostics.

%% Validate the one-run, three-parameter fit contract
    if ~isscalar(trace)
        error("TGS:OneRunFit", ...
            "fitTwoLayerTGS fits one run at a time; use fitTwoLayerTGSRuns for repeats.");
    end

    if isfield(S,"fitParameters")
        error("TGS:RemovedOption", ...
            "S.fitParameters has been removed; alpha_f, alpha_s, and R are always fitted.");
    end

    [x0,lowerX,upperX] = physicalSearchSpace(S);
    options = fittingOptions(S);

%% Fit a smoothed trace to initialize the physical parameters
    initializationTrace = smoothThermalTrace(trace,S);
    initializationWeight = traceWeight(initializationTrace);
    initializationData = initializationWeight*initializationTrace.y(:);
    dummyData = zeros(size(initializationData));
    initializationModel = @(x,xdata) weightedThermalModel( ...
        x,xdata,initializationTrace,S,initializationWeight);
    xThermal = lsqcurvefit(initializationModel,x0,dummyData, ...
        initializationData,lowerX,upperX,options);
    xThermal = xThermal(:).';

%% Fit the unsmoothed thermal and displacement response
    fitTrace = finalFitTrace(trace,S);
    finalStarts = unique([xThermal;x0],"rows","stable");
    finalWeight = traceWeight(fitTrace);
    measuredWeighted = finalWeight*fitTrace.y(:);
    dummyData = zeros(size(measuredWeighted));
    thermalModel = @(x,xdata) weightedThermalModel( ...
        x,xdata,fitTrace,S,finalWeight);

    startResults = repmat(struct( ...
        "x0",[],"x",[],"resnorm",Inf,"exitflag",NaN,"output",[]), ...
        size(finalStarts,1),1);
    bestResnorm = Inf;
    selectedStart = NaN;

    for startIndex = 1:size(finalStarts,1)
        [candidateX,candidateResnorm,candidateResidual, ...
            candidateExitflag,candidateOutput,~,candidateJacobian] = ...
            lsqcurvefit(thermalModel,finalStarts(startIndex,:), ...
            dummyData,measuredWeighted,lowerX,upperX,options);

        startResults(startIndex).x0 = finalStarts(startIndex,:);
        startResults(startIndex).x = candidateX(:).';
        startResults(startIndex).resnorm = candidateResnorm;
        startResults(startIndex).exitflag = candidateExitflag;
        startResults(startIndex).output = candidateOutput;

        if candidateResnorm < bestResnorm
            bestX = candidateX(:).';
            bestResnorm = candidateResnorm;
            objectiveResidual = candidateResidual;
            exitflag = candidateExitflag;
            output = candidateOutput;
            jacobian = candidateJacobian;
            selectedStart = startIndex;
        end
    end

    if ~isfinite(selectedStart)
        error("TGS:Optimization", ...
            "The optimizer did not return a finite thermal-model solution.");
    end

%% Reconstruct the unweighted thermal-model signal
    [yfit,rebuilt] = weightedThermalModel(bestX,[],fitTrace,S,1);
    measured = fitTrace.y(:);
    physicalValues = 10.^bestX;

%% Calculate uncertainty and identifiability diagnostics
    boundTolerance = S.dx;
    if isfield(S,"boundTolerance") && isfinite(S.boundTolerance) && ...
            S.boundTolerance >= 0
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
%PHYSICALSEARCHSPACE Validate and log-transform the three physical inputs.

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

function fitTrace = finalFitTrace(trace,S)
%FINALFITTRACE Apply an optional fixed post-pump delay to the final fit.

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
        error("TGS:FitWindow", ...
            "The final fit window contains fewer than 10 points.");
    end

    fitTrace = trace;
    fitTrace.t = t(use);
    fitTrace.y = y(use);
end

function options = fittingOptions(S)
%FITTINGOPTIONS Build bounded least-squares options for both thermal stages.

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
%SMOOTHTHERMALTRACE Smooth the full-resolution signal for initialization.

    if isfield(trace,"tFull") && isfield(trace,"yFull") && ...
            ~isempty(trace.tFull) && ~isempty(trace.yFull)
        sourceTime = trace.tFull(:);
        sourceSignal = trace.yFull(:);
    else
        sourceTime = trace.t(:);
        sourceSignal = trace.y(:);
    end

    timeSteps = diff(sourceTime);
    if numel(sourceTime) < 3 || any(~isfinite(sourceTime)) || ...
            any(~isfinite(sourceSignal)) || any(timeSteps <= 0)
        error("TGS:Trace", ...
            "The run must contain finite, increasing time samples.");
    end
    if ~isscalar(S.smoothTime) || ~isfinite(S.smoothTime) || ...
            S.smoothTime < 0
        error("TGS:Smoothing", ...
            "S.smoothTime must be a nonnegative finite duration.");
    end

    if S.smoothTime == 0
        smoothedSignal = sourceSignal;
    else
        smoothedSignal = movmean(sourceSignal,S.smoothTime, ...
            "SamplePoints",sourceTime,"Endpoints","shrink");
    end

    thermalTrace = trace;
    thermalTrace.t = trace.t(:);
    thermalTrace.y = interp1(sourceTime,smoothedSignal, ...
        thermalTrace.t,"linear");
end

function [weightedFit,state] = weightedThermalModel( ...
        x,~,trace,S,weight)
%WEIGHTEDTHERMALMODEL Profile temperature, displacement, and offset terms.

    p = physicalParameters(x);
    [temperature,displacement] = TwoLayerModel( ...
        trace.t,trace.Lambda,p,S);
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
%PHYSICALPARAMETERS Convert three log10 coordinates to physical values.

    x = x(:).';
    if numel(x) ~= 3 || any(~isfinite(x))
        error("TGS:Parameters", ...
            "The optimizer must supply log10([alpha_f, alpha_s, R]).");
    end
    p = 10.^x;
end

function [yfit,coefficients] = linearFit(design,y)
%LINEARFIT Solve scaled linear nuisance amplitudes at fixed physical values.

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
%TRACEWEIGHT Scale one residual by measured pre-pump noise when available.

    y = trace.y(:);
    if any(~isfinite(y))
        error("TGS:DataScale", ...
            "The measured run contains nonfinite values.");
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
        error("TGS:DataScale", ...
            "The measured run has no signal variation.");
    end
    weight = 1/scale;
end

function statistics = nonlinearStatistics(J,resnorm,linearParameterCount)
%NONLINEARSTATISTICS Calculate covariance and scaled-Jacobian rank metrics.

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
