function plot_run_vs_current_by_concentration(dataTable, analysisFolder, targetVgsStr)
    
    if ~exist(analysisFolder, 'dir')
        mkdir(analysisFolder);
    end

    targetVGS = str2double(targetVgsStr);

    if isnan(targetVGS)
        error('Invalid VGS value entered.');
        return;
    end

    runNames = unique(dataTable.RunName);
    runNames = sort_nat(runNames); % Sort run names numerically if applicable
    concentrations = unique(dataTable.Concentration);
    cmap = hsv(length(concentrations));

    runX = 1:length(runNames);
    concentrationData = struct();

    for c = 1:length(concentrations)
        conc = concentrations(c);
        avgCurrents = nan(1, length(runNames));
        stdCurrents = nan(1, length(runNames));

        for r = 1:length(runNames)
            runName = runNames{r};
            runData = dataTable(strcmp(dataTable.RunName, runName) & ...
                                dataTable.Concentration == conc, :);

            voltages = [];
            currents = [];

            for j = 1:height(runData)
                data = readmatrix(runData.FilePath{j});
                if size(data, 2) < 4, continue; end

                v = data(:, 3);
                i = data(:, 4);

                if isempty(voltages)
                    voltages = v(:)';
                    currents = i(:)';
                elseif length(v) == length(voltages)
                    currents(end+1, :) = i(:)';
                else
                    warning('Skipping mismatched file: %s', runData.FilePath{j});
                end
            end

            if isempty(currents), continue; end

            meanI = mean(currents, 1, 'omitnan');
            stdI = std(currents, 0, 1, 'omitnan');

            [~, idx] = min(abs(voltages - targetVGS));
            avgCurrents(r) = meanI(idx);
            stdCurrents(r) = stdI(idx);
        end

        concentrationData(c).Concentration = conc;
        concentrationData(c).Avg = avgCurrents;
        concentrationData(c).Std = stdCurrents;
    end

    % Combined Plot
    figure;
    hold on;
    legendEntries = {};
  
    for c = 1:length(concentrationData)
        errorbar(runX, concentrationData(c).Avg, concentrationData(c).Std, '-o', ...
                 'Color', cmap(c, :), 'LineWidth', 2);
        legendEntries{end+1} = sprintf('[%.2f]', concentrationData(c).Concentration);
    end

    xlabel('Run Number');
    ylabel(sprintf('Avg Current at V_{GS} = %.2f V', targetVGS));
    title('Run vs Current per Concentration');
    xticks(runX);
    xticklabels(runNames);
    legend(legendEntries, 'Location', 'eastoutside');
    grid on;
    hold off;

    saveas(gcf, fullfile(analysisFolder, 'Run_vs_Current_by_Concentration.png'));

    % Individual Plot and Excel per Concentration
    for c = 1:length(concentrationData)
        conc = concentrationData(c).Concentration;

        % Plot
        figure;
        errorbar(runX, concentrationData(c).Avg, concentrationData(c).Std, '-o', ...
                 'Color', cmap(c, :), 'LineWidth', 2);
        xlabel('Run Number');
        ylabel(sprintf('Avg Current at V_{GS} = %.2f V', targetVGS));
        title(sprintf('Run vs Current - [%.2f]', conc));
        xticks(runX);
        xticklabels(runNames);
        ylim([0.000 0.0035]);
        grid on;

        indivPlotFile = fullfile(analysisFolder, sprintf('Run_vs_Current_Conc_[%.2f].png', conc));
        saveas(gcf, indivPlotFile);

        % Excel for this concentration
        excelFile = fullfile(analysisFolder, sprintf('Run_vs_Current_Conc_[%.2f].xlsx', conc));
        headers = {'RunName', sprintf('[%.2f] Avg', conc), sprintf('[%.2f] Std', conc)};
        excelData = cell(length(runNames)+1, 3);
        excelData(1, :) = headers;
        for r = 1:length(runNames)
            excelData{r+1, 1} = runNames{r};
            excelData{r+1, 2} = concentrationData(c).Avg(r);
            excelData{r+1, 3} = concentrationData(c).Std(r);
        end
        writecell(excelData, excelFile);
    end

    % Save Combined Excel Summary
    summaryFile = fullfile(analysisFolder, 'Run_vs_Current_by_Concentration.xlsx');
    header = [{'RunName'}];
    for c = 1:length(concentrations)
        concStr = sprintf('[%.2f]', concentrations(c));
        header{end+1} = [concStr ' Avg'];
        header{end+1} = [concStr ' Std'];
    end

    outputCell = cell(length(runNames)+1, length(header));
    outputCell(1, :) = header;
    for r = 1:length(runNames)
        outputCell{r+1, 1} = runNames{r};
        col = 2;
        for c = 1:length(concentrationData)
            outputCell{r+1, col} = concentrationData(c).Avg(r);
            outputCell{r+1, col+1} = concentrationData(c).Std(r);
            col = col + 2;
        end
    end

    writecell(outputCell, summaryFile);
end

function sorted = sort_nat(c)
    % Natural-order sort for cell array of strings with numbers
    [~, idx] = sort(str2double(regexprep(c, '\D', '')));
    sorted = c(idx);
end

