function plot_perc_change_by_cell_name(dataTable, analysisFolder, targetVgsStr)
    % This function generates a plot for each unique cell, showing the
    % percentage change in current vs. concentration. The percentage change
    % is calculated relative to the average current at zero concentration.
    % It also generates a summary plot and table, averaging the results
    % for each letter group (e.g., 'A', 'B', 'C') from the CellName.
    %
    % Inputs:
    %   dataTable: A table containing columns 'CellName', 'Concentration',
    %              and 'FilePath'. It is assumed that there is at least one
    %              data point with a Concentration of 0 for each cell.
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
    
    % Initialize a cell array to store individual summary tables
    allSummaryTables = cell(length(cellNames), 1);
    
    % Loop through each unique cell name
    for c = 1:length(cellNames)
        cellName = cellNames{c};
        
        % Filter the data table to get only the data for the current cell
        cellData = dataTable(strcmp(dataTable.CellName, cellName), :);
        
        % --- Step 1: Find the baseline (zero concentration) data ---
        baseData = cellData(cellData.Concentration == 0, :);
        if isempty(baseData)
            warning('No concentration 0 data found for cell: %s. Skipping plot.', cellName);
            continue;
        end
        
        baseCurrentsAtVGS = [];
        % Read the zero concentration files to get the baseline current
        for j = 1:height(baseData)
            data = readmatrix(baseData.FilePath{j});
            if size(data, 2) < 4, continue; end
            v = data(:, 2);
            i = data(:, 4);
            [~, idx] = min(abs(v - targetVGS));
            baseCurrentsAtVGS(end+1) = i(idx);
        end
        
        % Calculate the average baseline current
        baseAvgCurrent = mean(baseCurrentsAtVGS, 'omitnan');
        if isnan(baseAvgCurrent) || baseAvgCurrent == 0
            warning('Baseline current for cell %s is zero or NaN. Skipping plot.', cellName);
            continue;
        end
        
        % --- Step 2: Process all concentrations and calculate percentage change ---
        
        % Get all unique concentrations for this cell, sorted
        concentrations = unique(cellData.Concentration);
        
        % Initialize arrays to store the summary data
        summaryConcs = [];
        summaryPercChange = [];
        summaryStdsPercChange = [];
        
        % Loop through each unique concentration
        for conc = concentrations'
            
            % Get all file paths for this specific concentration
            files = cellData(cellData.Concentration == conc, :).FilePath;
            
            runPercChanges = [];
            % Read data from all files for the current concentration
            for j = 1:height(files)
                data = readmatrix(files{j});
                if size(data, 2) < 4, continue; end
                v = data(:, 2);
                i = data(:, 4);
                % Find the current at the target VGS for this specific run
                [~, idx] = min(abs(v - targetVGS));
                currentAtVGS = i(idx);
                % Calculate the percentage change for this individual run
                percChange = ((currentAtVGS - baseAvgCurrent) / baseAvgCurrent) * 100;
                runPercChanges(end+1) = percChange;
            end
            
            % If no data was read, skip to the next concentration
            if isempty(runPercChanges), continue; end
            
            % Calculate the mean and standard deviation of the percentage changes
            meanPercChange = mean(runPercChanges, 'omitnan');
            stdPercChange = std(runPercChanges, 'omitnan');
            
            % Store the summary data
            summaryConcs(end+1) = conc;
            summaryPercChange(end+1) = meanPercChange;
            summaryStdsPercChange(end+1) = stdPercChange;
        end
        
        % Sort the summary data by concentration (if not already sorted)
        [summaryConcs, sortIdx] = sort(summaryConcs);
        summaryPercChange = summaryPercChange(sortIdx);
        summaryStdsPercChange = summaryStdsPercChange(sortIdx);
        
        % --- Step 3: Create plot and save output for the individual cell ---
        
        figure('Visible', 'off'); % Make the figure invisible to avoid pop-up windows
        errorbar(summaryConcs, summaryPercChange, summaryStdsPercChange, '-o', 'LineWidth', 2);
        xlabel('Concentration [pM]');
        ylabel('Current Percentage Change (%)');
      
        title(sprintf('Cell %s - Percentage Change vs Concentration', cellName));
        grid on;
        
        % Save the plot as a PNG file
        summaryPlotFile = fullfile(analysisFolder, sprintf('Cell_%s_perc_change_vs_concentration.png', cellName));
        saveas(gcf, summaryPlotFile);
        
        % Create the summary Excel file path
        summaryExcelFile = fullfile(analysisFolder, sprintf('Cell_%s_perc_change_vs_concentration.csv', cellName));
        
        % Create a summary table and save it as a CSV file
        summaryTable = table(summaryConcs(:), summaryPercChange(:), summaryStdsPercChange(:), ...
            'VariableNames', {'Concentration', 'AvgPercentageChange', 'StdDevPercentageChange'});
        writetable(summaryTable, summaryExcelFile);
        
        % Store the summary table in the cell array for later aggregation
        allSummaryTables{c} = summaryTable;
        
        % Close the current figure
        close(gcf);
    end

    % =========================================================================
    % %% New Section: Aggregate and Plot Data by Letter Group
    % This section processes the individual summary tables to create a single
    % plot and table that averages the data for all cells in the same letter
    % group (e.g., all "Pizza_A_Slice_X" cells are grouped together).
    % =========================================================================

    % Create a map to store all percentage change data for each letter group
    % Key: letter group (e.g., 'A')
    % Value: a struct with fields 'Concentrations' and 'PercChanges'
    groupData = containers.Map;
    
    % Loop through the cell names and their corresponding summary tables
    for c = 1:length(cellNames)
        cellName = cellNames{c};
        currentTable = allSummaryTables{c};
        
        % Skip if the table is empty (e.g., due to a warning in the main loop)
        if isempty(currentTable)
            continue;
        end
        
        % Extract the letter group from the cell name
        parts = strsplit(cellName, '_');
        if length(parts) > 1
            % Handle names with underscores (e.g., Pizza_A_Slice_X)
            group = parts{2}; 
        else
            % Handle names that are just a letter and a number (e.g., A1, B2)
            group = regexp(cellName, '^[A-Za-z]', 'match', 'once');
        end
        
        % If this is the first time we see this group, initialize it
        if ~isKey(groupData, group)
            groupData(group) = struct('Concentrations', currentTable.Concentration, ...
                                      'PercChanges', []);
        end
        
        % Append the percentage change data for the current cell
        currentGroup = groupData(group);
        % Reshape to a column vector for concatenation
        concatenatedData = [currentGroup.PercChanges, currentTable.AvgPercentageChange(:)];
        currentGroup.PercChanges = concatenatedData;
        groupData(group) = currentGroup;
    end
    
    % Get all unique letter groups
    letterGroups = keys(groupData);
    
    % Loop through each letter group and create a plot and summary table
    for k = 1:length(letterGroups)
        group = letterGroups{k};
        currentGroup = groupData(group);
        
        % Calculate the averaged percentage change and its standard deviation
        avgPercChange = mean(currentGroup.PercChanges, 2, 'omitnan');
        stdPercChange = std(currentGroup.PercChanges, 0, 2, 'omitnan');
        
        % --- Create plot for the averaged group data ---
        
        figure('Visible', 'off');
        errorbar(currentGroup.Concentrations, avgPercChange, stdPercChange, '-o', 'LineWidth', 2);
        xlabel('Concentration [pM]');
        ylabel('Average Current Percentage Change (%)');
        title(sprintf('Group %s - Averaged Percentage Change vs Concentration (n=%d)', group, size(currentGroup.PercChanges, 2)));
        grid on;
        
        % Save the plot as a PNG file
        avgPlotFile = fullfile(analysisFolder, sprintf('Group_%s_avg_perc_change_vs_concentration.png', group));
        saveas(gcf, avgPlotFile);
        
        % Create a summary table and save it as a CSV file
        avgSummaryTable = table(currentGroup.Concentrations(:), avgPercChange(:), stdPercChange(:), ...
            'VariableNames', {'Concentration', 'AvgPercentageChange', 'StdDevPercentageChange'});
        avgExcelFile = fullfile(analysisFolder, sprintf('Group_%s_avg_perc_change_vs_concentration.csv', group));
        writetable(avgSummaryTable, avgExcelFile);
        
        close(gcf);
    end
end