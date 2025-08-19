function plot_voltage_vs_current(dataTable, analysisFolder)

    if ~exist(analysisFolder, 'dir')
        mkdir(analysisFolder);
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
            cmap = hsv(length(concentrations));

            figure;
            hold on;
            legendEntries = {};
            avgData = [];
   
            for c = 1:length(concentrations)
                conc = concentrations(c);
                files = runData(runData.Concentration == conc, :).FilePath;

                voltages = [];
                currents = [];

                for j = 1:height(files)
                    data = readmatrix(files{j});
                    if size(data, 2) < 4, continue; end

                    % Modified: Change column 3 to 2 for x-axis data
                    v = data(:, 2);
                    i = data(:, 4);

                    if isempty(voltages)
                        voltages = v(:)';
                        currents = i(:)';
                    elseif length(v) == length(voltages)
                        currents(end+1, :) = i(:)';
                    else
                        warning('Skipping mismatched length file: %s', files{j});
                    end
                end

                if isempty(currents), continue; end

                meanI = mean(currents, 1, 'omitnan');
                stdI = std(currents, 0, 1, 'omitnan');

                fill([voltages, fliplr(voltages)], ...
                     [meanI + stdI, fliplr(meanI - stdI)], ...
                     cmap(c, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');

                h(c) = plot(voltages, meanI, 'Color', cmap(c, :), 'LineWidth', 2);
                legendEntries{end+1} = sprintf('[%.2f]', conc);

                avgData(c).allCurrents = currents;
                avgData(c).concentration = conc;
                avgData(c).voltage = voltages;
                avgData(c).meanCurrent = meanI;
                avgData(c).stdCurrent = stdI;
            end

            % Modified: Change x-axis label
            xlabel('Drain Voltage (V)');
            ylabel('Average Current [A]');
            title(sprintf('Cell: %s | %s', cellName, runName));
            legend(h, legendEntries, 'Location', 'eastoutside');
            set(gca, 'Position', [0.1, 0.1, 0.65, 0.8]);
            grid on;
            hold off;
           % ylim([0 0.00]);

            plotFile = fullfile(analysisFolder, sprintf('%s_%s_avg_plot.png', cellName, runName));
            saveas(gcf, plotFile);

            excelFile = fullfile(analysisFolder, sprintf('%s_%s_avg_data.xlsx', cellName, runName));
            dataMatrix = avgData(1).voltage(:);
            % Modified: Change VGS to VDS in the header
            headers = {'VDS'};

            for k = 1:length(avgData)
                conc = avgData(k).concentration;
                sweepCurrents = avgData(k).allCurrents;
                numSweeps = size(sweepCurrents, 1);

                for s = 1:numSweeps
                    dataMatrix(:, end+1) = sweepCurrents(s, :)';
                    headers{end+1} = sprintf('[%.2f] Sweep %d', conc, s);
                end

                dataMatrix(:, end+1) = avgData(k).meanCurrent(:);
                headers{end+1} = sprintf('[%.2f] Avg', conc);

                dataMatrix(:, end+1) = avgData(k).stdCurrent(:);
                headers{end+1} = sprintf('[%.2f] Std', conc);
            end

            outputCell = [headers; num2cell(dataMatrix)];
            writecell(outputCell, excelFile);
        end
    end
end