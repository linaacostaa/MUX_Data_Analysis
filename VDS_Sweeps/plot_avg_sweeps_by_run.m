% Author: Lina M. Acosta Perez
% Description: Creates a single plot for each measurement run, showing the average current vs. voltage 
% for all concentrations within that run. It also saves the summary data to an Excel file.
% Date: Sept. 12 2025

function plot_avg_sweeps_by_run(dataTable, masterFolder)
   
    if ~exist(masterFolder, 'dir')
        mkdir(masterFolder);
    end

    % Create a dedicated folder for plots and Excel files inside the master folder
    runSummaryFolder = fullfile(masterFolder, 'run_avg_sweeps');
    if ~exist(runSummaryFolder, 'dir')
        mkdir(runSummaryFolder);
    end

    runNames = unique(dataTable.RunName);

    for r = 1:length(runNames)
        runName = runNames{r};
        runData = dataTable(strcmp(dataTable.RunName, runName), :);

        figure;
        hold on;
        legendHandles = [];
        legendEntries = {};
      
        concentrations = unique(runData.Concentration);
        cmap = hsv(length(concentrations));

        % Modified: Change VGS to VDS for the new x-axis
        summaryTableData = {'VDS'};
        voltages = [];

        for c = 1:length(concentrations)
            conc = concentrations(c);
            concData = runData(runData.Concentration == conc, :);

            currents = [];

            for j = 1:height(concData)
                data = readmatrix(concData.FilePath{j});
                if size(data, 2) < 4
                    continue;
                end

                % Modified: Change column 3 to 2 for x-axis data
                v = data(:, 2);  % Drain Voltage
                i = data(:, 4);  % Current

                if isempty(voltages)
                    voltages = v(:)';
                end

                if length(v) == length(voltages)
                    currents(end+1, :) = i(:)';
                else
                    warning('Skipping mismatched length file: %s', concData.FilePath{j});
                end
            end

            if isempty(currents)
                continue;
            end

            meanI = mean(currents, 1, 'omitnan');
            stdI = std(currents, 0, 1, 'omitnan');

            % Plot with standard deviation shading
            fill([voltages, fliplr(voltages)], ...
                 [meanI + stdI, fliplr(meanI - stdI)], ...
                 cmap(c, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            hLine = plot(voltages, meanI, 'Color', cmap(c, :), 'LineWidth', 2);
            legendHandles(end+1) = hLine;
            legendEntries{end+1} = sprintf('[%.2f]', conc);

            % Append to summary table
            if c == 1
                summaryTableData(2:length(voltages)+1, 1) = num2cell(voltages');
            end
            summaryTableData{1, end+1} = sprintf('[%.2f] Avg Current [A]', conc);
            summaryTableData(2:length(meanI)+1, end) = num2cell(meanI');
            summaryTableData{1, end+1} = sprintf('[%.2f] Std Dev [A]', conc);
            summaryTableData(2:length(stdI)+1, end) = num2cell(stdI');
        end

        % Modified: Change x-axis label
        xlabel('Drain Voltage (V)');
        ylabel('Average Current [A]');
        title(sprintf('Average Current by Concentration - Run: %s', runName));
        legend(legendHandles, legendEntries, 'Location', 'eastoutside');
        set(gca, 'Position', [0.1, 0.1, 0.65, 0.8]);
        grid on;
        hold off;

        % Save the figure
        plotFile = fullfile(runSummaryFolder, sprintf('Run_%s_avg_sweeps.png', runName));
        saveas(gcf, plotFile);
        fprintf('✅ Saved average concentration plot for run %s: %s\n', runName, plotFile);

        % Save the Excel file
        excelFile = fullfile(runSummaryFolder, sprintf('Run_%s_avg_sweeps.xlsx', runName));
        writecell(summaryTableData, excelFile);
        fprintf('✅ Saved Excel for run %s: %s\n', runName, excelFile);
    end
end