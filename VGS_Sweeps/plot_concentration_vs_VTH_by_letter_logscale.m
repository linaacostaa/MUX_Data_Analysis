% Author: Lina M. Acosta Perez
% Description: Finds the average VGS at a specific targetCurrentStr for each concentration within a group (determined by the first letter of the CellName). It then plots concentration vs. average VGS with error bars. The x-axis is in logarithmic scale.
% Date: Sept. 12 2025

function plot_concentration_vs_VTH_by_letter_logscale(dataTable, analysisFolder, targetCurrentStr)
% plot_concentration_vs_current_by_group_letter_logscale Generates plots and saves data for concentration vs. VGS at a specific current, grouped by letter, with log scale on X-axis.
% This function prompts the user for a target current value. For each
% concentration within each group letter, it finds the average VGS at which
% the current is closest to the target current. It then plots Concentration
% vs. Average VGS with error bars (X-axis in log scale) and saves the plot
% as a PNG and the summary data as an Excel file.
%
% Inputs:
%   dataTable        - A MATLAB table containing the gathered data,
%                      expected to have columns like 'CellName', 'Concentration',
%                      and 'FilePath'. This function adds a 'GroupLetter' column.
%   analysisFolder   - The path to the folder where the generated plots
%                      and Excel files will be saved.

    % Ensure the analysis folder exists
    if ~exist(analysisFolder, 'dir')
        mkdir(analysisFolder);
    end

    targetCurrent = str2double(targetCurrentStr);
    
    fprintf('Processing for Target Current = %.2e A...\n', targetCurrent); % Informative message


    % Extract the first letter of each CellName to create a 'GroupLetter' column
    dataTable.GroupLetter = cellfun(@(s) regexp(s, '^[A-Za-z]', 'match', 'once'), ...
                                    dataTable.CellName, 'UniformOutput', false);

    % Get unique group letters from the data table
    groupLetters = unique(dataTable.GroupLetter);

    % Loop through each unique group letter
    for g = 1:length(groupLetters)
        group = groupLetters{g};
        % Filter data for the current group
        groupData = dataTable(strcmp(dataTable.GroupLetter, group), :);

        % Get unique concentrations for the current group
        concentrations = unique(groupData.Concentration);

        % Initialize arrays to store summary data for plotting
        summaryConcs = [];
        summaryVGSValues = []; % Stores average VGS at target current
        summaryVGSStds = [];   % Stores standard deviation of VGS at target current

        % Loop through each unique concentration within the group
        for c = 1:length(concentrations)
            conc = concentrations(c);
            % Get file paths for the current concentration and group
            files = groupData(groupData.Concentration == conc, :).FilePath;

            allVGSAtTargetCurrent = []; % To collect VGS values from each sweep at the target current

            % Loop through each file (sweep) associated with the current concentration and group
            for j = 1:height(files)
                % Read data from the file (assuming it's a numeric matrix)
                data = readmatrix(files{j});

                % Check if the file has enough columns (VGS is 3rd, Current is 4th)
                if size(data, 2) < 4
                    warning('Skipping file %s: Not enough columns (expected at least 4).', files{j});
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

            % If no valid VGS data was gathered for this concentration, skip to next
            if isempty(allVGSAtTargetCurrent)
                warning('No valid VGS data found for Concentration %f in Group %s at Target Current = %.2e A. Skipping.', conc, group, targetCurrent);
                continue;
            end

            % Calculate mean and standard deviation of VGS values found across sweeps
            meanVGS = mean(allVGSAtTargetCurrent, 'omitnan');
            stdVGS = std(allVGSAtTargetCurrent, 0, 'omitnan'); % 0 for default normalization

            % Store the concentration, mean VGS, and standard deviation of VGS
            summaryConcs(end+1) = conc;
            summaryVGSValues(end+1) = meanVGS;
            summaryVGSStds(end+1) = stdVGS;
        end

        % Sort concentrations and corresponding VGS values/stds for proper plotting
        [summaryConcs, sortIdx] = sort(summaryConcs);
        summaryVGSValues = summaryVGSValues(sortIdx);
        summaryVGSStds = summaryVGSStds(sortIdx);

        % Skip plotting if no summary data was collected
        if isempty(summaryConcs)
            warning('No summary data to plot for Group %s at Target Current = %.2e A. Skipping plot generation.', group, targetCurrent);
            continue;
        end

        % Create a new figure for each plot
        figure;
        % Plot concentration vs. VGS with error bars
        errorbar(summaryConcs, summaryVGSValues, summaryVGSStds, '-o', 'LineWidth', 2);

        % Add labels and title
        xlabel('Concentration [pM]');
        xscale log; % Maintain log scale on X-axis as per original function's name/intent
        ylabel(sprintf('Avg V_{GS} at Current = %.2e A', targetCurrent));
        title(sprintf('Group %s - Concentration vs V_{GS}', group));
        % ylim auto-scales now, removed fixed limit for VGS

        grid on; % Add a grid for better readability

        % Construct file names including the Current value to avoid overwrites
        current_str_raw = sprintf('Current_%.3e', targetCurrent);
      
        summaryPlotFile = fullfile(analysisFolder, sprintf('Group_%s_concentration_vs_VGS_%s.png', group, current_str_raw));
        summaryExcelFile = fullfile(analysisFolder, sprintf('Group_%s_concentration_vs_VGS_%s.xlsx', group, current_str_raw)); % Original was .xlsx

        % Save the plot as a PNG image
        saveas(gcf, summaryPlotFile);
        % Close the figure to prevent too many open figures
        close(gcf);

        % Create a table from summary data and save it to an Excel file
        summaryTable = table(summaryConcs(:), summaryVGSValues(:), summaryVGSStds(:), ...
                             'VariableNames', {'Concentration', 'AvgVGS', 'StdDevVGS'});
        writetable(summaryTable, summaryExcelFile);
    end
    fprintf('Finished plotting for the specified target current, grouped by letter (log scale).\n');
end
