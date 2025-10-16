% Author: Lina M. Acosta Perez
% Description: The purpose of these scripts is to automate the visualization and summarization of electrical data from multiple measurement runs
% and cells. The scripts generate various plots and summary tables to help users quickly assess device performance across different concentrations and runs.
% Date: Sept. 12 2025

clc; clear; close all;

%% Step 2: Ask the user to select the folder with concentration folders
parentFolderPath = uigetdir('', 'Select the parent folder that contains the [Concentration]/Run folders');
if parentFolderPath == 0
    error('No parent folder selected. Exiting...');
end

% Get all subfolders in the parent folder, which should be your concentration folders
subFolders = dir(parentFolderPath);
subFolders = subFolders([subFolders.isdir]);
subFolders = subFolders(~ismember({subFolders.name}, {'.', '..'}));

 %% Sampling Method Selection

    samplingChoice = questdlg('Which sampling method would you like to perform?', ...
        'Sampling Method', ...
        'IDS sampling', 'VDS sampling', 'IDS sampling');

    idsValue = '';
    vdsValue = '';

for i = 1:length(subFolders)
    concFolderPath = fullfile(parentFolderPath, subFolders(i).name);
    fprintf('Processing folder: %s\n', concFolderPath);
%% Gather data for concentration plots for the current folder
    dataTable = gather_data_for_concentration_plotting(concFolderPath);
 %% Step 4: Compile all CSV files into one massive table (massive compile)
% -------------------------------------------------------------------------
% This final version ensures the column headers maintain the 
% Parameter_Concentration_Run_FileName structure while strictly enforcing
% the length limit to prevent the array2table error.
% -------------------------------------------------------------------------

% Initialize variables
massiveCompileCell = {}; 
maxRows = 0;
allVarNames = {};

fprintf('\n--- Starting massive compile of CSV files (Side-by-Side Raw Data) ---\n');

% 1. Loop through all 'Parameter' folders directly under the parent
% NOTE: Assuming 'parentFolderPath' is defined earlier in your script.
paramFolders = dir(parentFolderPath); 
paramFolders = paramFolders([paramFolders.isdir] & ~ismember({paramFolders.name}, {'.', '..'}));

if isempty(paramFolders)
    warning('MATLAB:FolderSearch', 'No "Parameter" level folders found directly under the parent folder. Execution may fail if the folder structure is not standard.');
end

% Set the maximum length for the "Source Name" part of the header.
% The final header will be: SourceName + '_' + VariableName (e.g., '_VDS')
% We ensure the SourceName is short enough to allow space for the variable name suffix.
max_len_base = 60; 

% Explicitly loop through the known structure
for i = 1:length(paramFolders) % Changed from 'soapFolders'
    paramPath = fullfile(parentFolderPath, paramFolders(i).name); % Changed from 'soapPath'
    paramName = paramFolders(i).name; % Changed from 'soapName'
    
    % 2. Loop through 'Concentration' folders under the Parameter folder
    concFolders = dir(paramPath); 
    concFolders = concFolders([concFolders.isdir] & ~ismember({concFolders.name}, {'.', '..'}));
    
    for j = 1:length(concFolders)
        concPath = fullfile(paramPath, concFolders(j).name); 
        concName = concFolders(j).name;
        
        % 3. Loop through 'Run' folders under the Concentration folder
        runFolders = dir(concPath);
        runFolders = runFolders([runFolders.isdir] & ~ismember({runFolders.name}, {'.', '..'}));
        
        for l = 1:length(runFolders)
            runPath = fullfile(concPath, runFolders(l).name);
            runName = runFolders(l).name;
            
            % 4. Find all CSV files inside the Run folder
            currentCsvFiles = dir(fullfile(runPath, '*.csv'));
            
            % -----------------------------------------------------------------
            % CRITICAL FILTER: Remove hidden files, especially macOS metadata files (._*)
            currentCsvFiles = currentCsvFiles(~strncmp({currentCsvFiles.name}, '._', 2));
            % -----------------------------------------------------------------
            
            for k = 1:length(currentCsvFiles)
                filePath = fullfile(runPath, currentCsvFiles(k).name);
                fileName = currentCsvFiles(k).name;
                
                try
                    % Read the raw data
                    opts = detectImportOptions(filePath);
                    currentTable = readtable(filePath, opts);
                    currentRows = size(currentTable, 1);
                    
                    % --- Column Labeling and Concatenation ---
                    
                    % 1. Construct sourceName using folder names (Parameter, Conc, Run) and filename
                    baseSourceName = [paramName, '_', concName, '_', runName, '_', strrep(fileName, '.csv', '')]; % Changed from 'soapName'
                    
                    % 2. Sanitize the name: Replace illegal characters with '_'
                    baseSourceName = regexprep(baseSourceName, '[^a-zA-Z0-9_]', '_');
                    
                    % 3. Truncate the base name to the safe limit
                    if length(baseSourceName) > max_len_base
                        sourceName = baseSourceName(1:max_len_base);
                        fprintf('  - Warning: Source name truncated to %d chars.\n', max_len_base);
                    else
                        sourceName = baseSourceName;
                    end
                    
                    % Rename columns with the unique label
                    oldVarNames = currentTable.Properties.VariableNames;
                    newVarNames = cellfun(@(x) [sourceName, '_', x], oldVarNames, 'UniformOutput', false);
                    
                    % Store variable names for the header row
                    allVarNames = [allVarNames, newVarNames]; %#ok<AGROW>

                    % Convert table to cell array for robust concatenation
                    currentCell = table2cell(currentTable);
                    
                    if isempty(massiveCompileCell)
                        % First file starts the compile
                        massiveCompileCell = currentCell;
                        maxRows = currentRows;
                    else
                        % Handle row count mismatch (pad with NaNs)
                        if currentRows > maxRows
                            % Current table is longer: pad the existing massiveCompileCell
                            padding = cell(currentRows - maxRows, size(massiveCompileCell, 2));
                            padding(:) = {NaN};
                            massiveCompileCell = [massiveCompileCell; padding];
                            maxRows = currentRows;
                        elseif currentRows < maxRows
                            % Current table is shorter: pad the currentCell
                            padding = cell(maxRows - currentRows, size(currentCell, 2));
                            padding(:) = {NaN};
                            currentCell = [currentCell; padding];
                        end
                        
                        % Concatenate side-by-side (column-wise)
                        massiveCompileCell = [massiveCompileCell, currentCell];
                    end
                    
                    fprintf('  - Added: %s (Rows: %d)\n', sourceName, currentRows);
                    
                catch ME
                    warning('MATLAB:CSVReadError', 'Could not read CSV file %s. Error: %s', fileName, ME.message);
                end
            end % End of CSV file loop
        end % End of Run folder loop
    end % End of Concentration folder loop
end % End of Parameter folder loop


% --- Final Save Operation ---

outputFileName = 'massive_compile.csv';
outputFilePath = fullfile(parentFolderPath, outputFileName); 

if ~isempty(massiveCompileCell)
    % Verify the number of column names matches the compiled data
    numCols = size(massiveCompileCell, 2);
    if length(allVarNames) ~= numCols
        warning('MATLAB:VariableNamesMismatch', 'Column name count mismatch. Using default names for compiled columns.');
        allVarNames = cellfun(@(x) ['Var', num2str(x)], num2cell(1:numCols), 'UniformOutput', false);
    end
    
    % Prepend the column headers
    finalCell = [allVarNames; massiveCompileCell];
    
    % Use writecell for the most reliable saving of mixed-type data
    writecell(finalCell, outputFilePath);

    fprintf('\n✅ Successfully created **massive_compile.csv** with %d columns and %d rows, saved to: %s\n', size(massiveCompileCell, 2), size(massiveCompileCell, 1), outputFilePath);
else
    fprintf('\n⚠️ No CSV files were successfully compiled. massive_compile.csv was not created.\n');
end

fprintf('----------------------------------------------------------------\n');

%% Max value for i-v sweeps
% --- Add the following lines immediately after the massiveCompileCell is fully populated ---

% --- Find the Absolute Maximum Current Value by Column Name ---

if isempty(allVarNames)
    warning('MATLAB:MissingHeaders', 'The allVarNames variable is missing. Cannot identify current columns by name. Assuming default column 4, 8, 12, etc.');
    % Fallback to previous column indexing logic
    numCols = size(massiveCompileCell, 2);
    currentColumns = 4:4:numCols;
else
    % 1. Use the header names (allVarNames) to identify columns that end with "Id_A_"
    % `endsWith` creates a logical array (true for current columns, false otherwise).
    isCurrentColumn = endsWith(allVarNames, 'Id_A_');

    % 2. Convert the logical array to a vector of column indices (e.g., [4, 8, 12, ...])
    currentColumns = find(isCurrentColumn);
end

if ~isempty(currentColumns)
    % 3. Convert the numeric data in the identified current columns to a double matrix.
    % We convert only the necessary columns for efficiency.
    currentDataCell = massiveCompileCell(:, currentColumns);
    dataMatrix = cell2mat(currentDataCell);

    % 4. Linearize the matrix, take the absolute value, and find the single maximum value.
    maxAbsoluteCurrent = max(abs(dataMatrix(:)));
    minAbsoluteCurrent = min(abs(dataMatrix(:)));

    % Display the final result to the command window.
    fprintf('\n--- DATA SUMMARY ---\n');
    fprintf('Identified %d current columns (ending in "Id_A_").\n', length(currentColumns));
    fprintf('The **ABSOLUTE MAXIMUM CURRENT** value found across all datasets is: **%.4e A**\n', maxAbsoluteCurrent);
    fprintf('--------------------\n');
    fprintf('The **ABSOLUTE MINIMUM CURRENT** value found across all datasets is: **%.4e A**\n', minAbsoluteCurrent);
    fprintf('--------------------\n');
else
    % Handle the case where no current columns were found
    maxAbsoluteCurrent = 0; % Set to 0 to prevent errors in subsequent plotting functions
    minAbsoluteCurrent = 0;
    fprintf('\n⚠️ WARNING: No columns ending in "Id_A_" were found. Max current set to 0.\n');
end   
    %% Step 4: Create a folder to save concentration plots inside the selected folder
    analysisFolder_concentration = fullfile(concFolderPath, '4_voltage_vs_current');
    if ~exist(analysisFolder_concentration, 'dir')
        mkdir(analysisFolder_concentration);
    end
    %% Step 4.1: Generate average current plots by concentration - only averages the 5 sweeps of a cell name at one run
    plot_voltage_vs_current(dataTable, analysisFolder_concentration, maxAbsoluteCurrent, minAbsoluteCurrent);
    %% Step 7: Create a folder to save concentration plots inside the selected folder
    analysisFolder_volt_vs_avg_current_letter = fullfile(concFolderPath, '7_volt_vs_current_avg_letter');
    if ~exist(analysisFolder_volt_vs_avg_current_letter, 'dir')
        mkdir(analysisFolder_volt_vs_avg_current_letter);
    end
    %% Step 7.1: Generate average current plots by concentration
    plot_volt_vs_avg_current_by_group_letter(dataTable, analysisFolder_volt_vs_avg_current_letter, maxAbsoluteCurrent, minAbsoluteCurrent);
    %% Step 11: Create a folder to save concentration plots inside the selected folder
    analysisFolder_avg_concentration = fullfile(concFolderPath, '11_sweep_avg_per_run');
    if ~exist(analysisFolder_avg_concentration, 'dir')
        mkdir(analysisFolder_avg_concentration);
    end
    %% Step 11.1: Generate average of all the sweeps for a given run number
    plot_avg_sweeps_by_run(dataTable, analysisFolder_avg_concentration, maxAbsoluteCurrent, minAbsoluteCurrent);

    switch samplingChoice
        case 'VDS sampling'
            vdsValue = inputdlg('Enter a value for VDS sampling:', 'VDS Value');
            if isempty(vdsValue)
                error('VDS value not provided. Exiting...');
            end
            vdsValue = vdsValue{1}; % Extract value from cell array

            % --- Find IDS max at the requested VDS using massiveCompileCell / allVarNames ---
            vdsNum = str2double(vdsValue);
            if isnan(vdsNum)
                warning('Provided VDS ("%s") is not numeric. Using autoscaled IDS limits.', vdsValue);
                max_IDS_raw = maxAbsoluteCurrent;
                min_IDS_raw = minAbsoluteCurrent;
            else
                tol = 1e-6; % tolerance for matching VDS (adjust if needed)

                % identify VDS and IDS columns in the header list (try multiple common suffixes)
                vdsMask = endsWith(allVarNames, 'Vds_V_') | endsWith(allVarNames, 'Vds_V') | contains(allVarNames, 'Vds');
                idsMask = endsWith(allVarNames, 'Id_A') | endsWith(allVarNames, 'Id_A_') | endsWith(allVarNames, 'Id_A_') | contains(allVarNames, 'Id_A') | contains(allVarNames, 'Ids');

                vdsCols = find(vdsMask);
                idsCols = find(idsMask);

                matchedIdsValues = [];

                if isempty(vdsCols) || isempty(idsCols)
                    warning('Could not find VDS or IDS columns by header patterns. Falling back to maxAbsoluteCurrent.');
                    max_IDS_raw = maxAbsoluteCurrent;
                    min_IDS_raw = minAbsoluteCurrent;
                else
                    % scan each VDS column for rows that match the requested VDS
                    for vc = 1:length(vdsCols)
                        colIdx = vdsCols(vc);
                        % extract column as numeric, robust to numeric or string cells
                        colNum = nan(size(massiveCompileCell,1),1);
                        for rr = 1:size(massiveCompileCell,1)
                            val = massiveCompileCell{rr, colIdx};
                            if isnumeric(val)
                                colNum(rr) = double(val);
                            else
                                colNum(rr) = str2double(string(val));
                            end
                        end

                        matchRows = find(~isnan(colNum) & abs(colNum - vdsNum) <= tol);
                        if isempty(matchRows)
                            % allow small relative tolerance if none found
                            relTol = 1e-3;
                            matchRows = find(~isnan(colNum) & abs(colNum - vdsNum) <= max(relTol*abs(vdsNum), tol));
                        end

                        if ~isempty(matchRows)
                            % for every matching row, collect IDS values from all idsCols
                            for rr = matchRows(:)'
                                for ic = 1:length(idsCols)
                                    idValRaw = massiveCompileCell{rr, idsCols(ic)};
                                    if isnumeric(idValRaw)
                                        idNum = double(idValRaw);
                                    else
                                        idNum = str2double(string(idValRaw));
                                    end
                                    if ~isnan(idNum)
                                        matchedIdsValues(end+1) = idNum; %#ok<AGROW>
                                    end
                                end
                            end
                        end
                    end

                    if isempty(matchedIdsValues)
                        warning('No IDS values found at VDS = %.6g in compiled data. Falling back to maxAbsoluteCurrent.', vdsNum);
                        max_IDS_raw = maxAbsoluteCurrent;
                        min_IDS_raw = minAbsoluteCurrent;
                    else
                        % use absolute maximum as requested
                        max_IDS_raw = max(abs(matchedIdsValues(:)));
                        min_IDS_raw = min(abs(matchedIdsValues(:)));
                    end
                end
            end
            fprintf('Computed max_IDS_raw = %.4e A\n', max_IDS_raw);
            fprintf('Computed min_IDS_raw = %.4e A\n', min_IDS_raw);

            % VDS Sampling functions
            analysisFolder_concentration_vs_current = fullfile(concFolderPath, ['5_conc_vs_current_VDS_sampling_' vdsValue]);
            if ~exist(analysisFolder_concentration_vs_current, 'dir')
                mkdir(analysisFolder_concentration_vs_current);
            end
            plot_concentration_vs_current(dataTable, analysisFolder_concentration_vs_current, vdsValue, max_IDS_raw, min_IDS_raw);
            
            analysisFolder_concentration_vs_current_logscale = fullfile(concFolderPath, ['6_conc_vs_current_logscale_VDS_sampling_' vdsValue]);
            if ~exist(analysisFolder_concentration_vs_current_logscale, 'dir')
                mkdir(analysisFolder_concentration_vs_current_logscale);
            end
            plot_concentration_vs_current_logscale(dataTable, analysisFolder_concentration_vs_current_logscale, vdsValue, max_IDS_raw, min_IDS_raw);
        
            analysisFolder_conc_vs_avg_current_well = fullfile(concFolderPath, ['8_conc_vs_current_avg_letter_VDS_sampling_' vdsValue]);
            if ~exist(analysisFolder_conc_vs_avg_current_well, 'dir')
                mkdir(analysisFolder_conc_vs_avg_current_well);
            end
            plot_concentration_vs_current_by_group_letter(dataTable, analysisFolder_conc_vs_avg_current_well,vdsValue, max_IDS_raw, min_IDS_raw);
            
            analysisFolder_conc_vs_avg_current_well_logscale = fullfile(concFolderPath, ['9_conc_vs_current_avg_well_logscale_VDS_sampling_' vdsValue]);
            if ~exist(analysisFolder_conc_vs_avg_current_well_logscale, 'dir')
                mkdir(analysisFolder_conc_vs_avg_current_well_logscale);
            end
            plot_concentration_vs_current_by_group_letter_logscale(dataTable, analysisFolder_conc_vs_avg_current_well_logscale, vdsValue, max_IDS_raw, min_IDS_raw);
            
            analysisFolder_plot_run_vs_current_by_concentration = fullfile(concFolderPath, ['10_run_vs_IDS_VDS_sampling_' vdsValue]);
            if ~exist(analysisFolder_plot_run_vs_current_by_concentration, 'dir')
                mkdir(analysisFolder_plot_run_vs_current_by_concentration);
            end
            plot_run_vs_current_by_concentration(dataTable, analysisFolder_plot_run_vs_current_by_concentration, vdsValue, max_IDS_raw, min_IDS_raw);

            analysisFolder_concentration_vs_current_by_cell_name = fullfile(concFolderPath, ['20_conc_vs_current_per_cell_VDS_sampling_' vdsValue]);
            if ~exist(analysisFolder_concentration_vs_current_by_cell_name, 'dir')
                mkdir(analysisFolder_concentration_vs_current_by_cell_name);
            end
            plot_concentration_vs_current_by_cell_name(dataTable, analysisFolder_concentration_vs_current_by_cell_name, vdsValue, max_IDS_raw,min_IDS_raw);

            analysisFolder_perc_change_by_cell_name = fullfile(concFolderPath, ['21_conc_vs_current_per_change_V_sampling_' vdsValue]);
            if ~exist(analysisFolder_perc_change_by_cell_name, 'dir')
                mkdir(analysisFolder_perc_change_by_cell_name);
            end
            plot_perc_change_by_cell_name(dataTable, analysisFolder_perc_change_by_cell_name, vdsValue);
            
        otherwise
            disp('No sampling method selected or invalid choice. Exiting script.');
    end
end


%%FUNCTIONS

% Author: Lina M. Acosta Perez
% Description: Generates a separate I-V plot for each unique cell and run. It includes a 
% line for each concentration and saves a summary Excel file with all sweep and average data.
% Date: Sept. 12 2025

function plot_voltage_vs_current(dataTable, analysisFolder, maxAbsoluteCurrent, minAbsoluteCurrent)

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
            ylim([minAbsoluteCurrent*0.90 maxAbsoluteCurrent*1.1]);

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

% Author: Lina M. Acosta Perez
% Description: Groups the data by the first letter of the CellName (e.g., all 'A' cells together) and plots the average current vs. voltage for each concentration within that group.
% Date: Sept. 12 2025


function plot_volt_vs_avg_current_by_group_letter(dataTable, analysisFolder, maxAbsoluteCurrent, minAbsoluteCurrent)

    if ~exist(analysisFolder, 'dir')
        mkdir(analysisFolder);
    end

    % Extract the group letter from each CellName (e.g., A1 -> A)
    dataTable.GroupLetter = cellfun(@(s) regexp(s, '^[A-Za-z]', 'match', 'once'), ...
                                    dataTable.CellName, 'UniformOutput', false);
    groupLetters = unique(dataTable.GroupLetter);

    for g = 1:length(groupLetters)
        group = groupLetters{g};
        groupData = dataTable(strcmp(dataTable.GroupLetter, group), :);

        concentrations = unique(groupData.Concentration);
        cmap = hsv(length(concentrations));

        figure;
        hold on;
        legendEntries = {};
        rawExport = {};  % initialize cell array for CSV export
        % Modified: Change VGS to VDS in the header
        rawExport{1, 1} = 'VDS';
        colIdx = 2;

        % NEW: Initialize the cell array for the new summary CSV
        summaryExport = {};
        summaryExport{1, 1} = 'VDS';
        summaryColIdx = 2;

        for c = 1:length(concentrations)
            conc = concentrations(c);
            files = groupData(groupData.Concentration == conc, :).FilePath;

            voltages = [];
            currents = [];
            sweepLabels = {};

            for j = 1:height(files)
                data = readmatrix(files{j});
                if size(data, 2) < 4, continue; end

                % Modified: Change column 3 to 2 for x-axis data
                v = data(:, 2);  % Drain Voltage
                i = data(:, 4);  % Current

                if isempty(voltages)
                    voltages = v(:)';
                    currents = i(:)';
                elseif length(v) == length(voltages)
                    currents(end+1, :) = i(:)';
                else
                    warning('Skipping mismatched file: %s', files{j});
                end

                sweepLabels{end+1} = sprintf('[%.2f] Sweep %d', conc, j);
            end

            if isempty(currents), continue; end

            meanI = mean(currents, 1, 'omitnan');
            stdI = std(currents, 0, 1, 'omitnan');

            fill([voltages, fliplr(voltages)], ...
                 [meanI + stdI, fliplr(meanI - stdI)], ...
                 cmap(c, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');

            h(c) = plot(voltages, meanI, 'Color', cmap(c, :), 'LineWidth', 2);
            legendEntries{end+1} = sprintf('[%.2f]', conc);

            if c == 1
                rawExport(2:length(voltages)+1, 1) = num2cell(voltages');
                % NEW: Populate VDS column for the summary CSV
                summaryExport(2:length(voltages)+1, 1) = num2cell(voltages');
            end

            for s = 1:size(currents, 1)
                cellNames = groupData.CellName(groupData.Concentration == conc);
rawExport{1, colIdx} = sprintf('%s - %s', cellNames{s}, sweepLabels{s});
                rawExport(2:length(voltages)+1, colIdx) = num2cell(currents(s, :)');
                colIdx = colIdx + 1;
            end

            rawExport{1, colIdx} = sprintf('[%.2f] Avg', conc);
            rawExport{1, colIdx+1} = sprintf('[%.2f] Std', conc);
            rawExport(2:length(meanI)+1, colIdx) = num2cell(meanI');
            rawExport(2:length(stdI)+1, colIdx+1) = num2cell(stdI');
            colIdx = colIdx + 2;

            % NEW: Add average and standard deviation to the new summary CSV
            summaryExport{1, summaryColIdx} = sprintf('[%.2f] Avg', conc);
            summaryExport{1, summaryColIdx+1} = sprintf('[%.2f] Std', conc);
            summaryExport(2:length(meanI)+1, summaryColIdx) = num2cell(meanI');
            summaryExport(2:length(stdI)+1, summaryColIdx+1) = num2cell(stdI');
            summaryColIdx = summaryColIdx + 2;
        end

        % Modified: Change x-axis label
        xlabel('Drain Voltage (V)');
        ylabel('Average Current [A]');
        title(sprintf('Group %s - Sweep Plot', group));
        legend(h, legendEntries, 'Location', 'eastoutside');
        set(gca, 'Position', [0.1, 0.1, 0.65, 0.8]);
        grid on;
        hold off;
        ylim([minAbsoluteCurrent*0.90 maxAbsoluteCurrent*1.1]);

        % Save figure
        plotFile = fullfile(analysisFolder, sprintf('Group_%s_sweeps.png', group));
        saveas(gcf, plotFile);

        % Save raw data to CSV
        csvFile = fullfile(analysisFolder, sprintf('Group_%s_sweeps.csv', group));
        writecell(rawExport, csvFile);

        % NEW: Save the summary data to a separate CSV
        summaryCsvFile = fullfile(analysisFolder, sprintf('Group_%s_summary.csv', group));
        writecell(summaryExport, summaryCsvFile);
    end
end

% Author: Lina M. Acosta Perez
% Description: Creates a single plot for each measurement run, showing the average current vs. voltage 
% for all concentrations within that run. It also saves the summary data to an Excel file.
% Date: Sept. 12 2025

function plot_avg_sweeps_by_run(dataTable, masterFolder, maxAbsoluteCurrent, minAbsoluteCurrent)
   
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
        ylim([minAbsoluteCurrent*0.90 maxAbsoluteCurrent*1.1]);

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

function plot_concentration_vs_current(dataTable, analysisFolder,targetVgsStr, max_IDS_raw, min_IDS_raw)
      
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
            ylabel(sprintf('Avg Current at V_{GS} = %.2f V', targetVGS));
            title(sprintf('Cell: %s | %s', cellName, runName));
            ylim([min_IDS_raw*0.90 max_IDS_raw*1.1]);
         
            grid on;
            % Save plot and Excel
            
            % Construct file names including the VGS value to avoid overwrites
            vgs_str_raw = sprintf('VGS_%.2f', targetVGS);
            % Replace characters that might be problematic in filenames (e.g., '.', '+', '-')
            vgs_str = strrep(vgs_str_raw, '.', 'p');
            vgs_str = strrep(vgs_str, '+', '');
            vgs_str = strrep(vgs_str, '-', 'm');

            summaryPlotFile = fullfile(analysisFolder, sprintf('%s_%s_concentration_vs_current_%s.png', cellName, runName, vgs_str));
            saveas(gcf, summaryPlotFile);
            summaryExcelFile = fullfile(analysisFolder, sprintf('%s_%s_concentration_vs_current_%s.xlsx', cellName, runName, vgs_str));
            summaryTable = table(summaryConcs(:), summaryCurrents(:), summaryStds(:), 'VariableNames', {'Concentration', 'AvgCurrent', 'StdDev'});
            writetable(summaryTable, summaryExcelFile);
        end
    end

end

% Author: Lina M. Acosta Perez
% Description: Plots average current vs. concentration for each cell and run. The current is taken at a targetVgsStr (e.g., '0.5'). The x-axis (concentration) in a logarithmic scale.
% Date: Sept. 12 2025

function plot_concentration_vs_current_logscale(dataTable, analysisFolder, targetVgsStr, max_IDS_raw, min_IDS_raw)
   
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
            ylim([min_IDS_raw*0.90 max_IDS_raw*1.1]);
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

% Author: Lina M. Acosta Perez
% Description: Averages the data for all cells and runs that share the same group letter, then plots the average current vs. concentration with error bars.
% Date: Sept. 12 2025

function plot_concentration_vs_current_by_group_letter(dataTable, analysisFolder, targetVgsStr, max_IDS_raw, min_IDS_raw)
    
    if ~exist(analysisFolder, 'dir')
        mkdir(analysisFolder);
  end

    targetVGS = str2double(targetVgsStr);
    if isnan(targetVGS)
        error('Invalid VGS value entered.');
        return;
    end

    dataTable.GroupLetter = cellfun(@(s) regexp(s, '^[A-Za-z]', 'match', 'once'), ...
                                    dataTable.CellName, 'UniformOutput', false);
    groupLetters = unique(dataTable.GroupLetter);

    for g = 1:length(groupLetters)
        group = groupLetters{g};
        groupData = dataTable(strcmp(dataTable.GroupLetter, group), :);

        concentrations = unique(groupData.Concentration);
        summaryConcs = [];
        summaryCurrents = [];
        summaryStds = [];

        for c = 1:length(concentrations)
            conc = concentrations(c);
            files = groupData(groupData.Concentration == conc, :).FilePath;

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
        ylabel(sprintf('Avg Current at V_{GS} = %.2f V', targetVGS));
        ylim([min_IDS_raw*0.90 max_IDS_raw*1.1]);
      
        title(sprintf('Group %s - Concentration vs Current', group));
        grid on;

        summaryPlotFile = fullfile(analysisFolder, sprintf('Group_%s_concentration_vs_current.png', group));
        saveas(gcf, summaryPlotFile);

        summaryExcelFile = fullfile(analysisFolder, sprintf('Group_%s_concentration_vs_current.csv', group));
        summaryTable = table(summaryConcs(:), summaryCurrents(:), summaryStds(:), ...
            'VariableNames', {'Concentration', 'AvgCurrent', 'StdDev'});
        writetable(summaryTable, summaryExcelFile);
    end
end

% Author: Lina M. Acosta Perez
% Description: Averages the data for all cells and runs that share the same group letter, then plots the average current vs. concentration with error bars. The x-axis is in logarithmic scale.
% Date: Sept. 12 2025

function plot_concentration_vs_current_by_group_letter_logscale(dataTable, analysisFolder, targetVgsStr, max_IDS_raw, min_IDS_raw)
   
     if ~exist(analysisFolder, 'dir')
        mkdir(analysisFolder);
     end

    targetVGS = str2double(targetVgsStr);

    if isnan(targetVGS)
        error('Invalid VGS value entered.');
        return;
    end
    
    dataTable.GroupLetter = cellfun(@(s) regexp(s, '^[A-Za-z]', 'match', 'once'), ...
                                    dataTable.CellName, 'UniformOutput', false);
    groupLetters = unique(dataTable.GroupLetter);

    for g = 1:length(groupLetters)
        group = groupLetters{g};
        groupData = dataTable(strcmp(dataTable.GroupLetter, group), :);

        concentrations = unique(groupData.Concentration);
        summaryConcs = [];
        summaryCurrents = [];
        summaryStds = [];

        for c = 1:length(concentrations)
            conc = concentrations(c);
            files = groupData(groupData.Concentration == conc, :).FilePath;

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
        ylim([min_IDS_raw*0.90 max_IDS_raw*1.1]);
       
        title(sprintf('Group %s - Concentration vs Current', group));
        grid on;

        summaryPlotFile = fullfile(analysisFolder, sprintf('Group_%s_concentration_vs_current.png', group));
        saveas(gcf, summaryPlotFile);

        summaryExcelFile = fullfile(analysisFolder, sprintf('Group_%s_concentration_vs_current.xlsx', group));
        summaryTable = table(summaryConcs(:), summaryCurrents(:), summaryStds(:), ...
            'VariableNames', {'Concentration', 'AvgCurrent', 'StdDev'});
        writetable(summaryTable, summaryExcelFile);
    end
end

% Author: Lina M. Acosta Perez
% Description: Plots average current versus run number for each concentration at a specific target voltage (VGS). It generates both combined and individual plots, as well as Excel tables.
% Date: Sept. 12 2025

function plot_run_vs_current_by_concentration(dataTable, analysisFolder, targetVdsStr, max_IDS_raw, min_IDS_raw)
    
    if ~exist(analysisFolder, 'dir')
        mkdir(analysisFolder);
    end

    targetVDS = str2double(targetVdsStr);

    if isnan(targetVDS)
        error('Invalid VGS value entered.');
        return;
    end

    runNames = unique(dataTable.RunName);
    [~, idx] = sort(str2double(regexprep(runNames, '\\D', '')));
    runNames = runNames(idx);
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

                v = data(:, 2);
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

            [~, idx] = min(abs(voltages - targetVDS));
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

    xlabel('Run Name');
    ylabel(sprintf('Avg Current at V_{DS} = %.2f V', targetVDS));
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
        xlabel('Run Name');
        ylabel(sprintf('Avg Current at V_{DS} = %.2f V', targetVDS));
        title(sprintf('Run vs Current - [%.2f]', conc));
        xticks(runX);
        xticklabels(runNames);
        ylim([min_IDS_raw*0.90 max_IDS_raw*1.1]);
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

function plot_concentration_vs_current_by_cell_name(dataTable, analysisFolder, targetVgsStr, max_IDS_raw, min_IDS_raw)
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

                v = data(:, 2);
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
        ylim([min_IDS_raw*0.90 max_IDS_raw*1.10]); % Set a consistent y-axis limit for all plots
      
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
