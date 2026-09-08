function figureHandles = plotTwoLayerTGSResults(results,S)
% Create a separate fit-and-sensitivity figure for every run.
% Inputs:
%   results       - Per-run output from MainTwoLayerScript.
%   S             - Plot settings containing figureVisible.
% Output:
%   figureHandles - One figure handle per run, in acquisition order.

    fits = results.runFits;
    figureHandles = gobjects(numel(fits),1);

    for runIndex = 1:numel(fits)
        fit = fits(runIndex);
        trace = fit.traces;
        figureHandles(runIndex) = figure("Visible",S.figureVisible, ...
            "Name",sprintf("TGS spot %d, run %g", ...
            results.spot,trace.run),"NumberTitle","off");
        layout = tiledlayout(figureHandles(runIndex),2,1, ...
            "TileSpacing","compact","Padding","compact");
        title(layout,sprintf("TGS fit and sensitivity: spot %d, run %g", ...
            results.spot,trace.run))

        %% Plot this measured run and its selected multistart fit
        nexttile(layout)
        hold on
        plot(trace.t*1e6,trace.y,"-","Color",[0.45,0.45,0.45], ...
            "LineWidth",0.8,"DisplayName","Measured run")
        plot(trace.t*1e6,trace.yfit,"-","Color",[0,0.447,0.741], ...
            "LineWidth",1.2,"DisplayName","Thermal-model fit")
        xlabel("Pump-relative time (\mus)")
        ylabel("Detector signal")
        title(sprintf("%.4g \\mum, run %g", ...
            trace.Lambda*1e6,trace.run))
        legend("Location","best")
        grid on

        %% Plot this run's physical-parameter sensitivities
        nexttile(layout)
        parameterName = categorical(["alpha_f","alpha_s","R"]);
        sensitivity = max(fit.sensitivity(:),realmin);
        bar(parameterName,sensitivity)
        set(gca,"YScale","log")
        ylabel("Local sensitivity ||dr/dlog_{10}(p)||_2")
        grid on

        %% Apply consistent figure styling
        set(findall(figureHandles(runIndex),"-property","FontName"), ...
            "FontName","Cambria")
        set(findall(figureHandles(runIndex),"-property","FontSize"), ...
            "FontSize",10)
    end
end
