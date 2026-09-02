function [traces,selectedSpot] = prepareTGSData(S)
%PREPARETGSDATA Load and preprocess one spot's repeated TGS runs.
% Inputs:
%   S            - File, calibration, and preprocessing settings.
% Outputs:
%   traces       - One prepared trace per complete POS/NEG run.
%   selectedSpot - Spot identifier selected by the calibration file.

%% Read and validate the selected calibration file
    calListing = dir(S.calFile);

    if isempty(calListing)
        error("TGS:Calibration", "No post-processing calibration files were found.");
    end

    if numel(calListing) ~= 1
        error("TGS:Calibration", ...
            "Choose exactly one post-processing calibration file; found %d.", ...
            numel(calListing));
    end

    T = readtable(fullfile(calListing.folder,calListing.name), ...
        "FileType","text", ...
        "Delimiter"," ", ...
        "MultipleDelimsAsOne",true, ...
        "ReadVariableNames",true, ...
        "VariableNamingRule","preserve");
    runName = string(T.run_name);
    LambdaUm = calibrationSpacing(T);

    selectedSpot = calibrationSpot(runName,calListing.name);

%% Build the nominal-to-calibrated grating map
    nominalUm = nan(numel(runName),1);
    for i = 1:numel(runName)
        nominalUm(i) = calibrationNominal(runName(i));
    end

    cal = containers.Map("KeyType","char","ValueType","double");
    nominalValues = unique(nominalUm(isfinite(nominalUm)));
    for i = 1:numel(nominalValues)
        nominal = nominalValues(i);
        candidates = LambdaUm(nominalUm == nominal & ...
            isfinite(LambdaUm) & LambdaUm > 0);
        if isempty(candidates)
            continue
        end
        [~,closest] = min(abs(candidates-nominal));
        cal(key(nominal)) = candidates(closest)*1e-6;
    end

    calKeys = keys(cal);
    nominalList = str2double(string(calKeys));
    LambdaList = cellfun(@(mapKey) cal(mapKey),calKeys);
    [nominalList,order] = sort(nominalList);
    calibratedList = 1e6*LambdaList(order);
    disp(table(nominalList(:),calibratedList(:), ...
        'VariableNames',{'Nominal_um','Calibrated_um'}))

%% Discover and parse measurement files
    listing = dir(fullfile(S.dataDir,S.filePattern));
    fileTemplate = struct("nominal",[],"spot",[],"polarity","", ...
        "run",[],"baseline",false,"path","");
    files = repmat(fileTemplate,numel(listing),1);
    fileCount = 0;
    
    for i = 1:numel(listing)
        token = regexp(listing(i).name,S.fileRegex,"names","once");
        if isempty(token)
            continue
        end

        f = fileTemplate;
        f.nominal = str2double(token.nominal);
        f.polarity = string(token.polarity);
        f.run = str2double(token.run);
        f.baseline = contains(token.location,"baseline");
        if f.baseline
            f.spot = NaN;
        else
            f.spot = str2double(extractAfter(string(token.location),"spot"));
        end

        f.path = string(fullfile(listing(i).folder,listing(i).name));
        fileCount = fileCount+1;
        files(fileCount) = f;
    end
    files = files(1:fileCount);

%% Pair POS and NEG files for the calibration-selected spot
    measurement = find(~[files.baseline] & [files.spot] == selectedSpot);
    baseline = find([files.baseline]);
    if isempty(measurement)
        error("TGS:SpotSelection", ...
            "No measurement files for spot %d were found in the selected directory.", ...
            selectedSpot);
    end

    pairKey = strings(numel(measurement),1);
    for i = 1:numel(measurement)
        f = files(measurement(i));
        pairKey(i) = sprintf("%.12g|%d|%d",f.nominal,f.spot,f.run);
    end

    uniquePair = unique(pairKey,"stable");
    traceTemplate = struct("Lambda",[],"q",[],"t",[],"y",[], ...
        "tFull",[],"yFull",[],"t0",[],"nominal",[],"spot",[], ...
        "run",[],"noiseStd",[],"files",struct());
    traces = repmat(traceTemplate,numel(uniquePair),1);
    traceCount = 0;

%% Assemble each calibrated differential run
    for i = 1:numel(uniquePair)
        member = measurement(pairKey == uniquePair(i));
        ip = member([files(member).polarity] == "POS");
        in = member([files(member).polarity] == "NEG");

        if numel(ip) ~= 1 || numel(in) ~= 1
            continue
        end
    
        f = files(ip);
        mapKey = key(f.nominal);
        if ~isKey(cal,mapKey)
            continue
        end

        [tbp,ybp,bpPath] = baselineTrace( ...
            f.nominal,"POS",files,baseline,S);
        [tbn,ybn,bnPath] = baselineTrace( ...
            f.nominal,"NEG",files,baseline,S);

        [tp, yp] = readTrace(files(ip).path, S);
        [tn, yn] = readTrace(files(in).path, S);

        [tp,yp,ybp] = commonGrid(tp,yp,tbp,ybp);
        [tn,yn,ybn] = commonGrid(tn,yn,tbn,ybn);
        yp = yp-ybp;
        yn = yn-ybn;
        [t,yp,yn] = commonGrid(tp,yp,tn,yn);
        y = 0.5*(yp-yn);

        t0 = arrivalTime(t, y, S);
        t = t - t0;
        prePump = t < 0;

        if nnz(prePump) >= 3
            noiseStd = std(y(prePump),1);
        else
            noiseStd = NaN;
        end

        use = t >= S.tFit(1) & t <= S.tFit(2);
        g = traceTemplate;
        g.Lambda = cal(mapKey);
        g.q = 2*pi/g.Lambda;
        g.t = t(use);
        g.y = y(use);
        g.tFull = [];
        g.yFull = [];
        g.t0 = t0;
        g.nominal = f.nominal;
        g.spot = f.spot;
        g.run = f.run;
        g.noiseStd = noiseStd;

        g.files = struct("POS",files(ip).path,"NEG",files(in).path, ...
            "baselinePOS",bpPath,"baselineNEG",bnPath);

        traceCount = traceCount+1;
        traces(traceCount) = g;
    end

%% Preserve full traces and select nonlinear-fit samples
    traces = traces(1:traceCount);
    traces = finalizeRuns(traces,S);
end

function traces = finalizeRuns(traces,S)
%FINALIZERUNS Preserve full data and downsample only nonlinear-fit arrays.

    for i = 1:numel(traces)
        t = traces(i).t(:);
        y = traces(i).y(:);
        traces(i).tFull = t;
        traces(i).yFull = y;

        % Retain logarithmically distributed points for nonlinear fitting.
        if isfield(S,"nFitPoints") && isfinite(S.nFitPoints) && ...
                S.nFitPoints > 1 && numel(t) > S.nFitPoints
            use = unique(round(logspace(0,log10(numel(t)),S.nFitPoints)));
            traces(i).t = t(use);
            traces(i).y = y(use);
        else
            traces(i).t = t;
            traces(i).y = y;
        end
    end
end

function [t,y,path] = baselineTrace(nominal,polarity,files,baseline,S)
%BASELINETRACE Return an exact or spacing-interpolated baseline trace.
    records = baseline([files(baseline).polarity] == polarity);
    nominalList = unique([files(records).nominal]);
    exact = records([files(records).nominal] == nominal);

    if isscalar(exact)
        [t,y] = readTrace(files(exact).path,S);
        path = files(exact).path;
        return
    end

    lowerNominal = max(nominalList(nominalList < nominal));
    upperNominal = min(nominalList(nominalList > nominal));

    if isempty(lowerNominal) || isempty(upperNominal)
        error("TGS:Baseline", ...
            "Nominal spacing %.12g is not bracketed by baseline traces.", nominal);
    end

    lower = records([files(records).nominal] == lowerNominal);
    upper = records([files(records).nominal] == upperNominal);

    if numel(lower) ~= 1 || numel(upper) ~= 1
        error("TGS:Baseline", ...
            "Each baseline spacing must contain one %s trace.", polarity);
    end

    [tl,yl] = readTrace(files(lower).path,S);
    [tu,yu] = readTrace(files(upper).path,S);
    [t,yl,yu] = commonGrid(tl,yl,tu,yu);

    weight = (nominal - lowerNominal)/(upperNominal - lowerNominal);
    y = (1 - weight)*yl + weight*yu;
    path = [files(lower).path; files(upper).path];
end

function [t,y] = readTrace(path,S)
%READTRACE Read configured columns and discard rows containing nonfinite data.

    A = readmatrix(path,"NumHeaderLines",S.nHeader);
    t = A(:,S.tCol)*S.tScale;
    y = A(:,S.yCol);
    use = isfinite(t) & isfinite(y);
    t = t(use);
    y = y(use);
end

function [t,yp,yn] = commonGrid(tp,yp,tn,yn)
%COMMONGRID Interpolate any two traces onto their shared time interval.
    
    if isequal(tp, tn)
        t = tp;
        return
    end    
   
    % Determine the interval over which both polarities contain data
        t1 = max(min(tp), min(tn));
        t2 = min(max(tp), max(tn));        
    
    % Use the time vector with the larger median sample spacing
        if median(diff(tp)) >= median(diff(tn))
            t = tp(tp >= t1 & tp <= t2);
        else
            t = tn(tn >= t1 & tn <= t2);
        end
    
    % Interpolate both polarity signals onto the selected common grid
        yp = interp1(tp, yp, t, "linear");
        yn = interp1(tn, yn, t, "linear");
end

function t0 = arrivalTime(t,y,S)
%ARRIVALTIME Estimate pump arrival from curvature near the signal maximum.
    index = find(t >= S.tSearch(1) & t <= S.tSearch(2));
    ys = y(index);
    [~, imax] = max(ys);
    d2y = gradient(gradient(ys));
    [~, id2] = max(d2y(1:imax));
    i0 = max(index(1), index(id2) - S.nShift);
    t0 = t(i0);
end

function k = key(value)
%KEY Create a stable text key for a nominal grating value.

    k = sprintf("%.12g", value);
end

function spacingUm = calibrationSpacing(calibrationTable)
%CALIBRATIONSPACING Read calibrated spacing from supported table headers.

    variableName = string(calibrationTable.Properties.VariableNames);
    supportedName = ["grating_spacing_um","grating_value[um]"];
    selected = find(ismember(lower(variableName),lower(supportedName)),1);
    if isempty(selected)
        error("TGS:Calibration", ...
            "Calibration table must contain grating_spacing_um or grating_value[um].");
    end
    spacingUm = double(calibrationTable.(variableName(selected)));
end

function nominalUm = calibrationNominal(runName)
%CALIBRATIONNOMINAL Extract nominal micrometers before the spot identifier.

    token = regexp(char(runName), ...
        "(?<nominal>\d+(?:\.\d+)?)(?:um)?(?=[_-]spot\d+)", ...
        "names","once");
    if isempty(token)
        nominalUm = NaN;
    else
        nominalUm = str2double(token.nominal);
    end
end

function spot = calibrationSpot(runName,calibrationName)
%CALIBRATIONSPOT Require calibration rows to identify exactly one spot.

    spotByRow = nan(numel(runName),1);

    for i = 1:numel(runName)
        token = regexp(char(runName(i)), ...
            "(?:^|[_-])spot(?<spot>\d+)(?=[_-]|$)", ...
            "names","once");

        if ~isempty(token)
            spotByRow(i) = str2double(token.spot);
        end
    end

    if isempty(spotByRow) || any(~isfinite(spotByRow))
        error("TGS:CalibrationSpot", ...
            "Every calibration run_name must contain a spot identifier such as '_spot1'.");
    end

    spotList = unique(spotByRow);

    if numel(spotList) ~= 1
        error("TGS:CalibrationSpot", ...
            "Choose calibration input for exactly one spot; found spots: %s.", ...
            strjoin(string(spotList),", "));
    end

    spot = spotList;

    nameToken = regexp(calibrationName, ...
        "(?:^|[_-])spot(?<spot>\d+)(?=[_-]|$)","names","once");

    if ~isempty(nameToken) && str2double(nameToken.spot) ~= spot
        error("TGS:CalibrationSpot", ...
            "Calibration filename selects spot %s, but its run_name rows select spot %d.", ...
            nameToken.spot,spot);
    end
end
