function combined = fitTwoLayerTGSRuns(traces,S)
% Fit repeated runs independently and combine estimates.
% Inputs:
%   traces - Prepared runs from one spot and one grating.
%   S - Model and optimizer settings passed to fitTwoLayerTGS.
% Output:
%   combined - Spot-level estimates, per-run fits, and diagnostics.

%% Validate and order the repeated runs
    if isempty(traces)
        error("TGS:Data","At least one run is required for fitting.");
    end

    if isfield(S,"fitParameters")
        error("TGS:RemovedOption", "S.fitParameters has been removed; alpha_f, alpha_s, and R are always fitted.");
    end

    if numel(unique([traces.Lambda])) ~= 1
        error("TGS:OneGrating", "fitTwoLayerTGSRuns expects repeated runs from exactly one grating.");
    end

    % File-system enumeration is lexical; report runs in acquisition order.
    [~,runOrder] = sort([traces.run]);
    traces = traces(runOrder);
    runCount = numel(traces);

%% Fit each run independently
    fitCell = cell(runCount,1);
    for runIndex = 1:runCount
        fprintf("fitting run %g (%d of %d)\n", traces(runIndex).run,runIndex,runCount);
        fitCell{runIndex} = fitTwoLayerTGS(traces(runIndex),S);
    end

    fits = [fitCell{:}];
    parameterValue = vertcat(fits.p);
    parameterError = vertcat(fits.parameterError);
    parameterIdentifiable = vertcat(fits.parameterIdentifiable);

%% Combine the three physical parameters across valid runs
    combinedValue = nan(1,3);
    combinedError = nan(1,3);
    withinFitError = nan(1,3);
    betweenRunStd = nan(1,3);
    runCountUsed = zeros(1,3);

    for parameterIndex = 1:3
        valid = parameterIdentifiable(:,parameterIndex) & isfinite(parameterValue(:,parameterIndex));
        runCountUsed(parameterIndex) = nnz(valid);

        if any(valid)
            [combinedValue(parameterIndex),combinedError(parameterIndex), withinFitError(parameterIndex),betweenRunStd(parameterIndex)] = combineRunEstimates(parameterValue(valid,parameterIndex), parameterError(valid,parameterIndex));
        else
            % Keep a finite diagnostic mean, but leave the error undefined.
            available = isfinite(parameterValue(:,parameterIndex));
            if any(available)
                combinedValue(parameterIndex) = mean( parameterValue(available,parameterIndex));
                betweenRunStd(parameterIndex) = std( parameterValue(available,parameterIndex));
            end
        end
    end

%% Build the per-run summary table
    run = [traces.run].';
    exitflag = [fits.exitflag].';
    identifiable = [fits.identifiable].';
    runSummary = table(run,parameterValue(:,1),parameterError(:,1), parameterValue(:,2),parameterError(:,2), parameterValue(:,3),parameterError(:,3),exitflag,identifiable,'VariableNames',{'Run','AlphaF','AlphaFError','AlphaS', 'AlphaSError','R','RError','ExitFlag','Identifiable'});

%% Package combined estimates and diagnostics
    combined.fits = fits;
    combined.runSummary = runSummary;
    combined.runCount = runCount;
    combined.p = combinedValue;
    combined.x = log10(combinedValue);
    combined.parameterError = combinedError;
    combined.withinFitError = withinFitError;
    combined.betweenRunStd = betweenRunStd;
    combined.parameterRunCount = runCountUsed;
    combined.parameterIdentifiable = runCountUsed > 0;
    combined.identifiable = all(combined.parameterIdentifiable);
    combined.sensitivity = mean(vertcat(fits.sensitivity),1,"omitnan");
    combined.atLowerBound = any(vertcat(fits.atLowerBound),1);
    combined.atUpperBound = any(vertcat(fits.atUpperBound),1);
    combined.traces = [fits.traces];
    combined.r = vertcat(fits.r);
    combined.resnorm = sum(combined.r.^2);
    combined.objectiveResnorm = sum([fits.resnorm]);
    combined.exitflag = exitflag;
    combined.output = {fits.output};
    combined.startResults = {fits.startResults};
    combined.jacobianRank = [fits.jacobianRank];
    combined.jacobianCondition = [fits.jacobianCondition];
    combined.degreesOfFreedom = [fits.degreesOfFreedom];
    combined.physicalJacobianRank = [fits.physicalJacobianRank];
    combined.physicalJacobianCondition = [fits.physicalJacobianCondition];
    combined.physicalDegreesOfFreedom = [fits.physicalDegreesOfFreedom];
end

function [meanValue,totalError,withinError,betweenStd] = combineRunEstimates(values,errors)
% Apply the single-layer meanAndError convention

    values = values(:);
    errors = errors(:);
    runCount = numel(values);
    meanValue = mean(values);
    betweenStd = std(values);
    finiteError = isfinite(errors);

    if any(finiteError)
        withinError = sqrt(sum(errors(finiteError).^2))/runCount;
        totalError = hypot(withinError,betweenStd);
    else
        withinError = NaN;
        totalError = NaN;
    end
end
