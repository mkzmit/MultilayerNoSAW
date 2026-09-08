function runResults = fitTwoLayerTGSRuns(traces,S)
% Fit repeated runs independently without averaging their estimates.
% Inputs:
%   traces - Prepared runs from one spot and one grating.
%   S - Model and optimizer settings passed to fitTwoLayerTGS.
% Output:
%   runResults - Ordered per-run fits, values, and diagnostics.

%% Validate and order the repeated runs
    if isempty(traces)
        error("TGS:Data","At least one run is required for fitting.");
    end

    if isfield(S,"fitParameters")
        error("TGS:RemovedOption", ...
            "S.fitParameters has been removed; alpha_f, alpha_s, and R are always fitted.");
    end

    if numel(unique([traces.Lambda])) ~= 1
        error("TGS:OneGrating", ...
            "fitTwoLayerTGSRuns expects repeated runs from exactly one grating.");
    end

    % File-system enumeration is lexical; report runs in acquisition order.
    [~,runOrder] = sort([traces.run]);
    traces = traces(runOrder);
    runCount = numel(traces);

%% Fit each run independently
    fitCell = cell(runCount,1);
    for runIndex = 1:runCount
        fprintf("fitting run %g (%d of %d)\n", ...
            traces(runIndex).run,runIndex,runCount);
        fitCell{runIndex} = fitTwoLayerTGS(traces(runIndex),S);
    end

    fits = [fitCell{:}];
    run = [traces.run].';
    parameterValue = vertcat(fits.p);
    parameterError = vertcat(fits.parameterError);
    sensitivity = vertcat(fits.sensitivity);
    parameterIdentifiable = vertcat(fits.parameterIdentifiable);
    exitflag = [fits.exitflag].';
    identifiable = [fits.identifiable].';
    objectiveResnorm = [fits.resnorm].';
    selectedStart = [fits.selectedStart].';
    startCount = arrayfun(@(oneFit) numel(oneFit.startResults),fits).';
    successfulStartCount = [fits.successfulStartCount].';

%% Build the per-run summary table
    runSummary = table(run,parameterValue(:,1),parameterError(:,1), ...
        parameterValue(:,2),parameterError(:,2), ...
        parameterValue(:,3),parameterError(:,3), ...
        objectiveResnorm,exitflag,identifiable,selectedStart, ...
        startCount,successfulStartCount, ...
        'VariableNames',{'Run','AlphaF','AlphaFError','AlphaS', ...
        'AlphaSError','R','RError','ObjectiveResnorm','ExitFlag', ...
        'Identifiable','SelectedStart','StartsTried','SuccessfulStarts'});

%% Package only per-run values and diagnostics
    runResults.fits = fits;
    runResults.runSummary = runSummary;
    runResults.runCount = runCount;
    runResults.run = run;
    runResults.p = parameterValue;
    runResults.x = vertcat(fits.x);
    runResults.sensitivity = sensitivity;
    runResults.parameterError = parameterError;
    runResults.parameterIdentifiable = parameterIdentifiable;
    runResults.identifiable = identifiable;
    runResults.traces = [fits.traces];
    runResults.r = {fits.r}.';
    runResults.resnorm = cellfun(@(residual) sum(residual.^2),runResults.r);
    runResults.objectiveResidual = {fits.objectiveResidual}.';
    runResults.objectiveResnorm = objectiveResnorm;
    runResults.exitflag = exitflag;
    runResults.output = {fits.output}.';
    runResults.startResults = {fits.startResults}.';
    runResults.selectedStart = selectedStart;
    runResults.startCount = startCount;
    runResults.successfulStartCount = successfulStartCount;
    runResults.atLowerBound = vertcat(fits.atLowerBound);
    runResults.atUpperBound = vertcat(fits.atUpperBound);
    runResults.jacobianRank = [fits.jacobianRank].';
    runResults.jacobianCondition = [fits.jacobianCondition].';
    runResults.degreesOfFreedom = [fits.degreesOfFreedom].';
    runResults.physicalJacobianRank = [fits.physicalJacobianRank].';
    runResults.physicalJacobianCondition = ...
        [fits.physicalJacobianCondition].';
    runResults.physicalDegreesOfFreedom = ...
        [fits.physicalDegreesOfFreedom].';
end
