function results = fitTwoLayerTGSRuns(traces,S)
% Fit repeated runs independently without averaging their estimates.
% Inputs:
%   traces - Prepared runs from one spot and one grating.
%   S - Model and optimizer settings passed to fitTwoLayerTGS.
% Output:
%   results - Ordered table of fitted values and selected start indices.

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

%% Build per-run thermal initializers from one shared FFT peak estimate
    initializationTraces = fftThermalInitialization(traces,S);

%% Fit each run independently
    fitCell = cell(runCount,1);
    for runIndex = 1:runCount
        fitCell{runIndex} = fitTwoLayerTGS( ...
            traces(runIndex),S,initializationTraces(runIndex));
    end

    fits = [fitCell{:}];
    run = [traces.run].';
    parameterValue = vertcat(fits.p);
    selectedStart = [fits.selectedStart].';

%% Return only the requested per-run values
    results = table(run,parameterValue(:,1),parameterValue(:,2), ...
        parameterValue(:,3),selectedStart, ...
        'VariableNames',{'Run','AlphaF','AlphaS','R','SelectedStart'});
end
