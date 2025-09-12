% Author: Lina M. Acosta Perez
% Description: Plots average VGS vs. concentration for each cell and run at a specific targetCurrentStr (e.g., '1e-6'). The x-axis (concentration) in a logarithmic scale.
% Date: Sept. 12 2025

function plot_concentration_vs_VTH_logscale(dataTable, analysisFolder, targetCurrentStr)
% plot_concentration_vs_current Generates plots and saves data for concentration vs. VGS at a specific current.
% This function prompts the user for a target current value. For each
% concentration within each cell and run, it finds the average VGS at which
% the current is closest to the target current. It then plots Concentration
% vs. Average VGS with error bars and saves the plot as a PNG and the
% summary data as an Excel file.
%
% Inputs:
%   dataTable        - A MATLAB table containing the gathered data,
%                      expected to have columns like 'CellName', 'RunName',
%                      'Concentration', and 'FilePath'.
%   analysisFolder   - The path to the folder where the generated plots
%                      and Excel files will be saved.

    % Ensure the analysis folder exists
    if ~exist(analysisFolder, 'dir')
        mkdir(analysisFolder);
    end

     % Convert the input string to a numeric target current
    % Removed the check for isempty(idsValue) and idsValue{1}
    % as the master script now handles the input dialog and extraction
    targetCurrent = str2double(targetCurrentStr); 
    
    % Validate targetCurrent conversion
    if isnan(targetCurrent)
        error('Invalid target current value provided: %s. Please enter a numeric value.', targetCurrentStr);
    end

    fprintf('Processing for Target Current = %.2e A...\n', targetCurrent); % Informative message

    % Get unique cell names from the data table
    cellNames = unique(dataTable.CellName);

    % Loop through each unique cell name
    for i = 1:length(cellNames)
        cellName = cellNames{i};
        % Filter data for the current cell
        cellData = dataTable(strcmp(dataTable.CellName, cellName), :);

        % Get unique run names for the current cell
        runNames = unique(cellData.RunName);

        % Loop through each unique run name
        for r = 1:length(runNames)
            runName = runNames{r};
            % Filter data for the current run
            runData = cellData(strcmp(cellData.RunName, runName), :);

            % Get unique concentrations for the current run
            concentrations = unique(runData.Concentration);

            % Initialize arrays to store summary data for plotting
            summaryConcs = [];
            summaryVGSValues = []; % Stores average VGS at target current
            summaryVGSStds = [];   % Stores standard deviation of VGS at target current

            % Loop through each unique concentration
            for c = 1:length(concentrations)
                conc = concentrations(c);
                % Get file paths for the current concentration
                files = runData(runData.Concentration == conc, :).FilePath;

                allVGSAtTargetCurrent = []; % To collect VGS values from each sweep at the target current

                % Loop through each file (sweep) associated with the current concentration
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
                    warning('No valid VGS data found for Concentration %f in Cell %s, Run %s at Target Current = %.2e A. Skipping.', conc, cellName, runName, targetCurrent);
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
                warning('No summary data to plot for Cell %s, Run %s at Target Current = %.2e A. Skipping plot generation.', cellName, runName, targetCurrent);
                continue;
            end

            % Create a new figure for each plot
            figure;
            % Plot concentration vs. VGS with error bars
            errorbar(summaryConcs, summaryVGSValues, summaryVGSStds, '-o', 'LineWidth', 2);

            % Add labels and title
            xlabel('Concentration [pM]');
            xscale log; % Keep log scale on X-axis as per your original code
            ylabel(sprintf('Avg V_{GS} at Current = %.2e A', targetCurrent));
            title(sprintf('Cell: %s | %s', cellName, runName));
            % ylim auto-scales now, removed fixed limit for VGS

            grid on; % Add a grid for better readability

            % Construct file names including the Current value to avoid overwrites
            current_str_raw = sprintf('Current_%.3e', targetCurrent);


            summaryPlotFile = fullfile(analysisFolder, sprintf('%s_%s_concentration_vs_VGS_%s.png', cellName, runName, current_str_raw));
            summaryExcelFile = fullfile(analysisFolder, sprintf('%s_%s_concentration_vs_VGS_%s.xlsx', cellName, runName, current_str_raw));

            % Save the plot as a PNG image
            saveas(gcf, summaryPlotFile);
          

            % Create a table from summary data and save it to an Excel file
            summaryTable = table(summaryConcs(:), summaryVGSValues(:), summaryVGSStds(:), ...
                                 'VariableNames', {'Concentration', 'AvgVGS', 'StdDevVGS'});
            writetable(summaryTable, summaryExcelFile);
        end
    end
    fprintf('Finished plotting for the specified target current.\n');
end
