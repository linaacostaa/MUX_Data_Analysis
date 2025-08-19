function plot_concentration_vs_current_by_cell_name(dataTable, analysisFolder, targetVgsStr)
    % This function generates an error bar plot for each unique cell,
    % showing the average current vs. concentration. The average and
    % standard deviation are calculated from all runs for a given
    % cell and concentration.
    %
    % Inputs:
    %   dataTable: A table containing columns 'CellName', 'Concentration',
    %              and 'FilePath'.
    %   analysisFolder: The directory to save the plots and summary data.
    %   targetVgsStr: A string representing the target gate-source voltage (VGS)
    %                 at which to extract the current data.
    
    % Create the analysis folder if it doesn't exist
    if ~exist(analysisFolder, 'dir')
        mkdir(analysisFolder);
    end

    % Convert the target VGS string to a numeric value
    targetVGS = str2double(targetVgsStr);
    if isnan(targetVGS)
        error('Invalid VGS value entered. It must be a number.');
        return;
    end

    % Find all unique cell names in the data table
    cellNames = unique(dataTable.CellName);

    % Loop through each unique cell name
    for c = 1:length(cellNames)
        cellName = cellNames{c};
        
        % Filter the data table to get only the data for the current cell
        cellData = dataTable(strcmp(dataTable.CellName, cellName), :);

        % Get the unique concentrations for this cell
        concentrations = unique(cellData.Concentration);
        
        % Initialize arrays to store the summary data for this cell
        summaryConcs = [];
        summaryCurrents = [];
        summaryStds = [];

        % Loop through each unique concentration for this cell
        for conc = concentrations'
            
            % Get all file paths for this specific cell and concentration
            files = cellData(cellData.Concentration == conc, :).FilePath;

            voltages = [];
            currents = [];

            % Read data from all files for the current concentration
            for j = 1:height(files)
                data = readmatrix(files{j});
                if size(data, 2) < 4, continue; end

                v = data(:, 3);
                i = data(:, 4);

                % Initialize the voltages and currents array with the first file
                if isempty(voltages)
                    voltages = v(:)';
                    currents = i(:)';
                % If subsequent files have the same voltage sweep length, add the current data
                elseif length(v) == length(voltages)
                    currents(end+1, :) = i(:)';
                else
                    warning('Skipping mismatched file: %s', files{j});
                end
            end

            % If no data was read, skip to the next concentration
            if isempty(currents), continue; end

            % Calculate the mean and standard deviation of the current across all runs
            meanI = mean(currents, 1, 'omitnan');
            stdI = std(currents, 0, 1, 'omitnan');

            % Find the index corresponding to the target VGS
            [~, idx] = min(abs(voltages - targetVGS));
            
            % Store the summary data for this concentration
            summaryConcs(end+1) = conc;
            summaryCurrents(end+1) = meanI(idx);
            summaryStds(end+1) = stdI(idx);
        end

        % Sort the summary data by concentration
        [summaryConcs, sortIdx] = sort(summaryConcs);
        summaryCurrents = summaryCurrents(sortIdx);
        summaryStds = summaryStds(sortIdx);

        % Create a new figure and plot the data with error bars
        figure;
        errorbar(summaryConcs, summaryCurrents, summaryStds, '-o', 'LineWidth', 2);
        xlabel('Concentration [pM]');
        ylabel(sprintf('Avg Current at V_{GS} = %.2f V', targetVGS));
        ylim([0 0.0035]); % Set a consistent y-axis limit for all plots
      
        title(sprintf('Cell %s - Concentration vs Current', cellName));
        grid on;

        % Save the plot as a PNG file
        summaryPlotFile = fullfile(analysisFolder, sprintf('Cell_%s_concentration_vs_current.png', cellName));
        saveas(gcf, summaryPlotFile);

        % Create a summary table and save it as a CSV file
        summaryExcelFile = fullfile(analysisFolder, sprintf('Cell_%s_concentration_vs_current.csv', cellName));
        summaryTable = table(summaryConcs(:), summaryCurrents(:), summaryStds(:), ...
            'VariableNames', {'Concentration', 'AvgCurrent', 'StdDev'});
        writetable(summaryTable, summaryExcelFile);
        
        % Close the current figure to prevent it from cluttering the display
        close(gcf);
    end
end
