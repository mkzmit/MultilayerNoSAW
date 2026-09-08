function figureHandle = plotTwoLayerTGSResults(results,S)
% Plot repeated-run fits and physical sensitivities.
% Inputs:
%   results      - Spot-level output from MainTwoLayerScript.
%   S            - Plot settings containing figureVisible.
% Output:
%   figureHandle - Handle to the generated figure.

%% Create the one-spot, one-grating layout
    traces = results.traces;
    figureHandle = figure("Visible",S.figureVisible,"Name",sprintf("TGS spot %d",results.spot));
    layout = tiledlayout(2,1,"TileSpacing","compact","Padding","compact");
    title(layout,sprintf("TGS fit and sensitivity for spot %d",results.spot))

%% Plot every measured run and its thermal-model fit
    nexttile
    hold on
    for runIndex = 1:numel(traces)
        trace = traces(runIndex);
        if runIndex == 1
            visibility = "on";
        else
            visibility = "off";
        end

        plot(trace.t*1e6,trace.y,"-","Color",[0.72,0.72,0.72], ...
            "LineWidth",0.7,"DisplayName","Measured runs", ...
            "HandleVisibility",visibility)
        plot(trace.t*1e6,trace.yfit,"-","Color",[0,0.447,0.741], ...
            "LineWidth",1,"DisplayName","Thermal-model fits", ...
            "HandleVisibility",visibility)
    end

    xlabel("Pump-relative time (\mus)")
    ylabel("Detector signal")
    title(sprintf("%.4g \\mum, %d runs", ...
        traces(1).Lambda*1e6,numel(traces)))
    legend("Location","best")
    grid on

%% Plot sensitivity for physical parameters
    nexttile
    parameterName = categorical(["alpha_f","alpha_s","R"]);
    sensitivity = max(results.sensitivity(:),realmin);
    bar(parameterName,sensitivity)
    set(gca,"YScale","log")
    ylabel("Local sensitivity ||dr/dlog_{10}(p)||_2")
    grid on

%% Apply consistent figure styling
    set(findall(figureHandle,"-property","FontName"), ...
        "FontName","Cambria")
    set(findall(figureHandle,"-property","FontSize"),"FontSize",10)
end
