% Author: Lina M. Acosta Perez
% Description: Plots average current vs. concentration for each cell and run. The current is taken at a targetVgsStr (e.g., '0.5'). The x-axis (concentration) in a logarithmic scale.
% Date: Sept. 12 2025

function plot_concentration_vs_current_logscale(dataTable, analysisFolder, targetVgsStr)
   
   if ~exist(analysisFolder, 'dir')
        mkdir(analysisFolder);
   end

    targetVGS = str2double(targetVgsStr);

    if isnan(targetVGS)
        error('Invalid VGS value entered.');
        return;
    end
    
    cellNames = unique(dataTable.CellName);

    for i = 1:length(cellNames)
        cellName = cellNames{i};
        cellData = dataTable(strcmp(dataTable.CellName, cellName), :);
        runNames = unique(cellData.RunName);

        for r = 1:length(runNames)
            runName = runNames{r};
            runData = cellData(strcmp(cellData.RunName, runName), :);

            concentrations = unique(runData.Concentration);
            summaryConcs = [];
            summaryCurrents = [];
            summaryStds = [];

            for c = 1:length(concentrations)
                conc = concentrations(c);
                files = runData(runData.Concentration == conc, :).FilePath;

                voltages = [];
                currents = [];

                for j = 1:height(files)
                    data = readmatrix(files{j});
                    if size(data, 2) < 4, continue; end

                    v = data(:, 2);
                    i = data(:, 4);

                    if isempty(voltages)
                        voltages = v(:)';
                        currents = i(:)';
                    elseif length(v) == length(voltages)
                        currents(end+1, :) = i(:)';
                    else
                        warning('Skipping mismatched file: %s', files{j});
                    end
                end

                if isempty(currents), continue; end

                meanI = mean(currents, 1, 'omitnan');
                stdI = std(currents, 0, 1, 'omitnan');

                [~, idx] = min(abs(voltages - targetVGS));
                summaryConcs(end+1) = conc;
                summaryCurrents(end+1) = meanI(idx);
                summaryStds(end+1) = stdI(idx);
            end

            [summaryConcs, sortIdx] = sort(summaryConcs);
            summaryCurrents = summaryCurrents(sortIdx);
            summaryStds = summaryStds(sortIdx);

            figure;
            errorbar(summaryConcs, summaryCurrents, summaryStds, '-o', 'LineWidth', 2);
            xlabel('Concentration [pM]');
            xscale log;
            ylabel(sprintf('Avg Current at V_{GS} = %.2f V', targetVGS));
            title(sprintf('Cell: %s | %s', cellName, runName));
            ylim([0 0.05]);
            grid on;
          
            % Save plot and Excel
            summaryPlotFile = fullfile(analysisFolder, sprintf('%s_%s_concentration_vs_current.png', cellName, runName));
            saveas(gcf, summaryPlotFile);

            summaryExcelFile = fullfile(analysisFolder, sprintf('%s_%s_concentration_vs_current.xlsx', cellName, runName));
            summaryTable = table(summaryConcs(:), summaryCurrents(:), summaryStds(:), 'VariableNames', {'Concentration', 'AvgCurrent', 'StdDev'});
            writetable(summaryTable, summaryExcelFile);
        end
    end
end