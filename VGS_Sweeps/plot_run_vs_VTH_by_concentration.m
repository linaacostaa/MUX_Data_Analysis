function plot_run_vs_VTH_by_concentration(dataTable, analysisFolder, targetCurrentStr)
% plot_run_vs_current_by_concentration Generates plots and saves data for Run Number vs. VGS at a specific current, grouped by concentration.
% This function prompts the user for a target current value. For each
% concentration, it finds the average VGS at which the current is closest
% to the target current across different runs. It then plots Run Number vs.
% Average VGS with error bars and saves combined and individual plots as
% PNGs, and summary data as Excel files.
%
% Inputs:
%   dataTable        - A MATLAB table containing the gathered data,
%                      expected to have columns like 'RunName', 'Concentration',
%                      and 'FilePath'.
%   analysisFolder   - The path to the folder where the generated plots
%                      and Excel files will be saved.

    if ~exist(analysisFolder, 'dir')
        mkdir(analysisFolder);
    end

    targetCurrent = str2double(targetCurrentStr);
    
    fprintf('Processing for Target Current = %.2e A...\n', targetCurrent); % Informative message

    runNames = unique(dataTable.RunName);
    runNames = sort_nat(runNames); % Sort run names numerically if applicable
    concentrations = unique(dataTable.Concentration);
    cmap = hsv(length(concentrations)); % Colormap for different concentrations

    runX = 1:length(runNames); % Numeric representation of run numbers for plotting

    concentrationData = struct(); % Structure to hold data for each concentration

    % Loop through each unique concentration
    for c = 1:length(concentrations)
        conc = concentrations(c);
        avgVGS = nan(1, length(runNames)); % Stores average VGS for each run at this concentration
        stdVGS = nan(1, length(runNames)); % Stores std dev of VGS for each run at this concentration

        % Loop through each run name
        for r = 1:length(runNames)
            runName = runNames{r};
            % Filter data for the current run and concentration
            runData = dataTable(strcmp(dataTable.RunName, runName) & ...
                                dataTable.Concentration == conc, :);

            allVGSAtTargetCurrent = []; % To collect VGS values from each sweep at the target current

            % Loop through each file (sweep) for the current run and concentration
            for j = 1:height(runData)
                data = readmatrix(runData.FilePath{j});

                % Check if the file has enough columns (VGS is 3rd, Current is 4th)
                if size(data, 2) < 4
                    warning('Skipping file %s: Not enough columns (expected at least 4).', runData.FilePath{j});
                    continue;
                end

                v = data(:, 3); % VGS (Gate Voltage)
                i = data(:, 4); % Current

                % Find the index where the current is closest to the target current
                [~, current_idx] = min(abs(i - targetCurrent));

                % Get the VGS value at that closest current point
                vgs_at_this_current = v(current_idx);

                % Store this VGS value
                allVGSAtTargetCurrent(end+1) = vgs_at_this_current;
            end

            % If valid VGS data was gathered for this run and concentration
            if ~isempty(allVGSAtTargetCurrent)
                avgVGS(r) = mean(allVGSAtTargetCurrent, 'omitnan');
                stdVGS(r) = std(allVGSAtTargetCurrent, 0, 'omitnan'); % 0 for default normalization
            else
                warning('No valid VGS data found for Run %s, Concentration %f at Target Current = %.2e A. Skipping.', runName, conc, targetCurrent);
            end
        end
        concentrationData(c).Concentration = conc;
        concentrationData(c).Avg = avgVGS;
        concentrationData(c).Std = stdVGS;
    end

    % Combined Plot for all concentrations
    figure;
    hold on;
    legendEntries = {};

    for c = 1:length(concentrationData)
        % Only plot if there is valid data for this concentration
        if ~all(isnan(concentrationData(c).Avg))
            errorbar(runX, concentrationData(c).Avg, concentrationData(c).Std, '-o', ...
                     'Color', cmap(c, :), 'LineWidth', 2);
            legendEntries{end+1} = sprintf('[%.2f]', concentrationData(c).Concentration);
        end
    end
    xlabel('Run Number');
    ylabel(sprintf('Avg V_{GS} at Current = %.2e A', targetCurrent));
    title('Run vs V_{GS} per Concentration');
    xticks(runX);
    xticklabels(runNames);
    legend(legendEntries, 'Location', 'eastoutside');
    grid on;
    hold off;

    % Construct file name including the Current value
    current_str_raw_combined = sprintf('Current_%.3e', targetCurrent);
   
    saveas(gcf, fullfile(analysisFolder, sprintf('Run_vs_VGS_by_Concentration_%s.png', current_str_raw_combined)));
    close(gcf); % Close figure after saving

    % Individual Plot and Excel per Concentration
    for c = 1:length(concentrationData)
        conc = concentrationData(c).Concentration;

        % Skip if no valid data for this concentration
        if all(isnan(concentrationData(c).Avg))
            warning('No valid data to generate individual plot/Excel for Concentration %.2f at Target Current = %.2e A. Skipping.', conc, targetCurrent);
            continue;
        end

        % Plot
        figure;
        errorbar(runX, concentrationData(c).Avg, concentrationData(c).Std, '-o', ...
                 'Color', cmap(c, :), 'LineWidth', 2);
        xlabel('Run Number');
        ylabel(sprintf('Avg V_{GS} at Current = %.2e A', targetCurrent));
        title(sprintf('Run vs V_{GS} - [%.2f]', conc));
        xticks(runX);
        xticklabels(runNames);
        % Removed fixed ylim as VGS range can vary
        grid on;

        % Construct file name including the Current value
        current_str_raw_indiv = sprintf('Current_%.3e', targetCurrent);


        indivPlotFile = fullfile(analysisFolder, sprintf('Run_vs_VGS_Conc_[%.2f]_%s.png', conc, current_str_raw_indiv));
        saveas(gcf, indivPlotFile);
        close(gcf); % Close figure after saving

        % Excel for this concentration
        excelFile = fullfile(analysisFolder, sprintf('Run_vs_VGS_Conc_[%.2f]_%s.xlsx', conc, current_str_raw_indiv));
        headers = {'RunName', sprintf('[%.2f] Avg VGS', conc), sprintf('[%.2f] Std VGS', conc)};
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
    summaryFile = fullfile(analysisFolder, sprintf('Run_vs_VGS_by_Concentration_%s.xlsx', current_str_raw_combined));
    header = [{'RunName'}];
    for c = 1:length(concentrations)
        concStr = sprintf('[%.2f]', concentrations(c));
        header{end+1} = [concStr ' Avg VGS'];
        header{end+1} = [concStr ' Std VGS'];
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
    writetable(cell2table(outputCell(2:end,:), 'VariableNames', outputCell(1,:)), summaryFile, 'WriteVariableNames', true);
    fprintf('Finished plotting for the specified target current (Run vs VGS by Concentration).\n');
end

function sorted = sort_nat(c)
    % Natural-order sort for cell array of strings with numbers
    % This helper function remains unchanged.
    [~, idx] = sort(str2double(regexprep(c, '\D', '')));
    sorted = c(idx);
end
