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
        'IDS sampling', 'VGS sampling', 'IDS sampling');

    idsValue = '';
    vgsValue = '';

for i = 1:length(subFolders)
    concFolderPath = fullfile(parentFolderPath, subFolders(i).name);
    fprintf('Processing folder: %s\n', concFolderPath);
%% Step 3: Gather data for concentration plots
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
% The final header will be: SourceName + '_' + VariableName (e.g., '_VGS')
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
    dataVector = dataMatrix(:);

    % Max absolute current (used for the upper limit)
    maxAbsoluteCurrent = max(abs(dataMatrix(:)));

    % Raw minimum current (used for the lower limit)
    minAbsoluteCurrent = min(dataVector); 

    % Display the final result to the command window.
    fprintf('\n--- DATA SUMMARY ---\n');
    fprintf('Identified %d current columns (ending in "Id_A_").\n', length(currentColumns));
    fprintf('The **ABSOLUTE MAXIMUM CURRENT** value found across all datasets is: **%.4e A**\n', maxAbsoluteCurrent);
    fprintf('--------------------\n');
    fprintf('The **RAW MINIMUM CURRENT** value found across all datasets is: **%.4e A**\n', minAbsoluteCurrent);
    fprintf('--------------------\n');
else
    % Handle the case where no current columns were found
    maxAbsoluteCurrent = 0; % Set to 0 to prevent errors in subsequent plotting functions
    minAbsoluteCurrent = 0;
    fprintf('\n⚠️ WARNING: No columns ending in "Id_A_" were found. Max/Min current set to 0.\n');
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
    case 'IDS sampling'
        % 1. Prompt user for IDS value
        idsValue = inputdlg('Enter a value for IDS sampling:', 'IDS Value');
        if isempty(idsValue)
            error('IDS value not provided. Exiting...');
        end
        idsValue = idsValue{1}; % Extract value from cell array

        % --- Setup and Conversion ---
        idsNum = str2double(string(idsValue)); 
        maxAbsoluteVoltage = 5; % Fallback max voltage
        tol = 1e-6; % Absolute tolerance for matching IDS

        if isnan(idsNum)
            warning('Provided IDS ("%s") is not numeric. Using autoscaled VGS limits.', idsValue);
            max_VGS_raw = maxAbsoluteVoltage;
            min_VGS_raw = 0;
        else
            % --- Stricter Column Identification (Kept from original code) ---
            vgsMask = endsWith(allVarNames, 'Vgs_V_') | endsWith(allVarNames, 'Vgs_V');
            idsMask = endsWith(allVarNames, 'Id_A_') | endsWith(allVarNames, 'Id_A');

            vgsCols = find(vgsMask);
            idsCols = find(idsMask);

            if isempty(vgsCols) || isempty(idsCols)
                warning('Could not find VGS or IDS columns using strict header patterns. Falling back to maxAbsoluteVoltage.');
                max_VGS_raw = maxAbsoluteVoltage;
                min_VGS_raw = 0;
            else
                % Data rows start from row 2 (assuming row 1 is the header)
                dataRows = size(massiveCompileCell, 1) - 1; 
                
                % Extract ONLY the data (excluding header row 1)
                massiveData = massiveCompileCell(2:end, :); 

                % --- FIX 1: Convert ALL relevant columns to numeric arrays (Vectorized) ---
                
                % Initialize numeric arrays for VGS and IDS
                VGS_numeric = nan(dataRows, length(vgsCols));
                IDS_numeric = nan(dataRows, length(idsCols));

                % Convert VGS columns
                for vc = 1:length(vgsCols)
                    colIdx = vgsCols(vc);
                    % Use cellfun to convert the entire column from cell to double
                    VGS_numeric(:, vc) = cellfun(@(x) str2double(string(x)), massiveData(:, colIdx));
                end

                % Convert IDS columns
                for ic = 1:length(idsCols)
                    colIdx = idsCols(ic);
                    IDS_numeric(:, ic) = cellfun(@(x) str2double(string(x)), massiveData(:, colIdx));
                end
                
                % --- FIX 2: Vectorized IDS Matching ---
                
                % Define a dynamic tolerance
                relTol = 1e-3;
                currentTol = max(relTol * abs(idsNum), tol);
                
                allMatchRows = [];
                
                % Iterate through each IDS column to find matching rows
                for ic = 1:length(idsCols)
                    current_ids_col = IDS_numeric(:, ic);
                    
                    % Vectorized finding of rows matching idsNum within tolerance
                    matchMask = ~isnan(current_ids_col) & (abs(current_ids_col - idsNum) <= currentTol);
                    
                    % Get row indices (1 to dataRows)
                    matchRows = find(matchMask);
                    allMatchRows = unique([allMatchRows; matchRows]); % Combine unique match rows
                end

                if isempty(allMatchRows)
                    warning('No VGS values found at IDS = %.6g in compiled data within tolerance.', idsNum);
                    max_VGS_raw = maxAbsoluteVoltage; % Fallback 
                else
                    % Extract all VGS values for the matched rows across ALL VGS columns
                    matchedVgsData = VGS_numeric(allMatchRows, :);
                    
                    % Flatten the array and remove NaN values
                    matchedVgsValues = matchedVgsData(~isnan(matchedVgsData));
                    
                    if isempty(matchedVgsValues)
                        warning('VGS values were NaN for all matched rows at IDS = %.6g.', idsNum);
                        max_VGS_raw = maxAbsoluteVoltage; 
                    else
                        % Calculate the maximum absolute VGS value
                        max_VGS_raw = max(abs(matchedVgsValues(:)));
                        min_VGS_raw = min(matchedVgsValues(:)); % New output
                    end
                end
            end
        
        fprintf('Computed max_Vgs_raw = %.4e V\n', max_VGS_raw);
        fprintf('Computed min_Vgs_raw = %.4e V\n', min_VGS_raw);
        % IDS Sampling functions
        analysisFolder_concentration_vs_VTH = fullfile(concFolderPath, ['12_conc_vs_VTH_IDS_sampling_' idsValue]);
        if ~exist(analysisFolder_concentration_vs_VTH, 'dir')
            mkdir(analysisFolder_concentration_vs_VTH);
        end
        plot_concentration_vs_VTH(dataTable, analysisFolder_concentration_vs_VTH, idsValue, max_VGS_raw, min_VGS_raw);
        
        analysisFolder_concentration_vs_VTH_logscale = fullfile(concFolderPath, ['13_conc_vs_VTH_logscale_IDS_sampling_' idsValue]);
        if ~exist(analysisFolder_concentration_vs_VTH_logscale, 'dir')
            mkdir(analysisFolder_concentration_vs_VTH_logscale);
        end
        plot_concentration_vs_VTH_logscale(dataTable, analysisFolder_concentration_vs_VTH_logscale, idsValue, max_VGS_raw, min_VGS_raw);
        
        analysisFolder_conc_vs_avg_VTH_well = fullfile(concFolderPath, ['14_conc_vs_VTH_by_letter_IDS_sampling_' idsValue]);
        if ~exist(analysisFolder_conc_vs_avg_VTH_well, 'dir')
            mkdir(analysisFolder_conc_vs_avg_VTH_well);
        end
        plot_concentration_vs_VTH_by_letter(dataTable, analysisFolder_conc_vs_avg_VTH_well, idsValue, max_VGS_raw, min_VGS_raw);
        
        analysisFolder_conc_vs_avg_VTH_well_logscale = fullfile(concFolderPath, ['15_conc_vs_VTH_by_letter_log_IDS_sampling_' idsValue]);
        if ~exist(analysisFolder_conc_vs_avg_VTH_well_logscale, 'dir')
            mkdir(analysisFolder_conc_vs_avg_VTH_well_logscale);
        end
        plot_concentration_vs_VTH_by_letter_logscale(dataTable, analysisFolder_conc_vs_avg_VTH_well_logscale, idsValue, max_VGS_raw, min_VGS_raw);
        
        analysisFolder_plot_run_vs_VTH_by_concentration = fullfile(concFolderPath, ['16_run_vs_VTH_IDS_sampling_' idsValue]);
        if ~exist(analysisFolder_plot_run_vs_VTH_by_concentration, 'dir')
            mkdir(analysisFolder_plot_run_vs_VTH_by_concentration);
        end
        plot_run_vs_VTH_by_concentration(dataTable, analysisFolder_plot_run_vs_VTH_by_concentration, idsValue, max_VGS_raw, min_VGS_raw);
    end
    %% VGS sampling
    case 'VGS sampling'
        vgsValue = inputdlg('Enter a value for VGS sampling:', 'VGS Value');
        if isempty(vgsValue)
            error('VGS value not provided. Exiting...');
        end
        vgsValue = vgsValue{1}; % Extract value from cell array

        % --- Find IDS max at the requested VGS using massiveCompileCell / allVarNames ---
        vgsNum = str2double(vgsValue);
        if isnan(vgsNum)
            warning('Provided VGS ("%s") is not numeric. Using autoscaled IDS limits.', vgsValue);
            max_IDS_raw = maxAbsoluteCurrent;
        else
            tol = 1e-6; % tolerance for matching VGS (adjust if needed)

            % identify VGS and IDS columns in the header list (try multiple common suffixes)
            vgsMask = endsWith(allVarNames, 'Vgs_V_') | endsWith(allVarNames, 'Vgs_V') | contains(allVarNames, 'Vgs');
            idsMask = endsWith(allVarNames, 'Id_A') | endsWith(allVarNames, 'Id_A_') | endsWith(allVarNames, 'Id_A_') | contains(allVarNames, 'Id_A') | contains(allVarNames, 'Ids');

            vgsCols = find(vgsMask);
            idsCols = find(idsMask);

            matchedIdsValues = [];

            if isempty(vgsCols) || isempty(idsCols)
                warning('Could not find VGS or IDS columns by header patterns. Falling back to maxAbsoluteCurrent.');
                max_IDS_raw = maxAbsoluteCurrent;
            else
                % scan each VGS column for rows that match the requested VGS
                for vc = 1:length(vgsCols)
                    colIdx = vgsCols(vc);
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

                    matchRows = find(~isnan(colNum) & abs(colNum - vgsNum) <= tol);
                    if isempty(matchRows)
                        % allow small relative tolerance if none found
                        relTol = 1e-3;
                        matchRows = find(~isnan(colNum) & abs(colNum - vgsNum) <= max(relTol*abs(vgsNum), tol));
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
                    warning('No IDS values found at VGS = %.6g in compiled data. Falling back to maxAbsoluteCurrent.', vgsNum);
                    max_IDS_raw = maxAbsoluteCurrent;
                    min_IDS_raw = 0;
                else
                    % use absolute maximum as requested
                    max_IDS_raw = max(abs(matchedIdsValues(:)));
                    min_IDS_raw = min(matchedIdsValues(:));
                end
            end
        end
        fprintf('Computed max_IDS_raw = %.4e A\n', max_IDS_raw);
        fprintf('Computed min_IDS_raw = %.4e A\n', min_IDS_raw); % New output

        % VGS Sampling functions
        analysisFolder_concentration_vs_current = fullfile(concFolderPath, ['5_conc_vs_current_VGS_sampling_' vgsValue]);
        if ~exist(analysisFolder_concentration_vs_current, 'dir')
            mkdir(analysisFolder_concentration_vs_current);
        end
        plot_concentration_vs_current(dataTable, analysisFolder_concentration_vs_current, vgsValue, max_IDS_raw, min_IDS_raw);
        
        analysisFolder_concentration_vs_current_logscale = fullfile(concFolderPath, ['6_conc_vs_current_logscale_VGS_sampling_' vgsValue]);
        if ~exist(analysisFolder_concentration_vs_current_logscale, 'dir')
            mkdir(analysisFolder_concentration_vs_current_logscale);
        end
        plot_concentration_vs_current_logscale(dataTable, analysisFolder_concentration_vs_current_logscale, vgsValue, max_IDS_raw, min_IDS_raw);
    
        analysisFolder_conc_vs_avg_current_well = fullfile(concFolderPath, ['8_conc_vs_current_avg_letter_VGS_sampling_' vgsValue]);
        if ~exist(analysisFolder_conc_vs_avg_current_well, 'dir')
            mkdir(analysisFolder_conc_vs_avg_current_well);
        end
        plot_concentration_vs_current_by_group_letter(dataTable, analysisFolder_conc_vs_avg_current_well,vgsValue, max_IDS_raw, min_IDS_raw);
        
        analysisFolder_conc_vs_avg_current_well_logscale = fullfile(concFolderPath, ['9_conc_vs_current_avg_well_logscale_VGS_sampling_' vgsValue]);
        if ~exist(analysisFolder_conc_vs_avg_current_well_logscale, 'dir')
            mkdir(analysisFolder_conc_vs_avg_current_well_logscale);
        end
        plot_concentration_vs_current_by_group_letter_logscale(dataTable, analysisFolder_conc_vs_avg_current_well_logscale, vgsValue, max_IDS_raw, min_IDS_raw);
        
        analysisFolder_plot_run_vs_current_by_concentration = fullfile(concFolderPath, ['10_parameter_vs_IDS_VGS_sampling_' vgsValue]);
        if ~exist(analysisFolder_plot_run_vs_current_by_concentration, 'dir')
            mkdir(analysisFolder_plot_run_vs_current_by_concentration);
        end
        plot_run_vs_current_by_concentration(dataTable, analysisFolder_plot_run_vs_current_by_concentration, vgsValue, max_IDS_raw, min_IDS_raw);

        analysisFolder_concentration_vs_current_by_cell_name = fullfile(concFolderPath, ['20_conc_vs_current_per_cell_VGS_sampling_' vgsValue]);
        if ~exist(analysisFolder_concentration_vs_current_by_cell_name, 'dir')
            mkdir(analysisFolder_concentration_vs_current_by_cell_name);
        end
        plot_concentration_vs_current_by_cell_name(dataTable, analysisFolder_concentration_vs_current_by_cell_name, vgsValue, max_IDS_raw, min_IDS_raw);

        analysisFolder_perc_change_by_cell_name = fullfile(concFolderPath, ['21_conc_vs_current_per_change_VGS_sampling_' vgsValue]);
        if ~exist(analysisFolder_perc_change_by_cell_name, 'dir')
            mkdir(analysisFolder_perc_change_by_cell_name);
        end
        plot_perc_change_by_cell_name(dataTable, analysisFolder_perc_change_by_cell_name, vgsValue);
        
    otherwise
        disp('No sampling method selected or invalid choice. Exiting script.');
end
end





%% Functions
% Author: Lina M. Acosta Perez
% Description: This function is responsible for navigating the data folder structure and compiling the master 
% dataTable used for all subsequent analysis. It automates the process of finding raw data files and their corresponding metadata.
% Date: Sept. 12 2025

function dataTable = gather_data_for_concentration_plotting(mainFolderPath)
% GATHER_DATA_FOR_CONCENTRATION_PLOTTING Collects data from a
% hierarchical folder structure for plotting.
% It now uses the new file naming convention instead of a JSON file.

dataRows = {};  % Initialize storage for table data

% Get concentration folders (e.g., [10], [50])
concentrationFolders = dir(mainFolderPath);
concentrationFolders = concentrationFolders([concentrationFolders.isdir]);
concentrationFolders = concentrationFolders(~ismember({concentrationFolders.name}, {'.', '..'}));

% Filter out folders starting with '.'
concentrationFolders = concentrationFolders(arrayfun(@(x) x.name(1) ~= '.', concentrationFolders));

for c = 1:length(concentrationFolders)
    concName = concentrationFolders(c).name;
    concPath = fullfile(mainFolderPath, concName);
    
    % Extract numeric concentration value from folder name [X]
    concMatch = regexp(concName, '\[(\d+\.?\d*)\]', 'tokens');
    if isempty(concMatch), continue; end
    concentration = str2double(concMatch{1}{1});
    
    % Get run folders inside this concentration folder (e.g., 1, 2, 3)
    runFolders = dir(concPath);
    runFolders = runFolders([runFolders.isdir]);
    runFolders = runFolders(~ismember({runFolders.name}, {'.', '..'}));
    
    % Filter out folders starting with '.'
    runFolders = runFolders(arrayfun(@(x) x.name(1) ~= '.', runFolders));
    
    for r = 1:length(runFolders)
        runName = runFolders(r).name;
        runPath = fullfile(concPath, runName);
        
        % Find data files in the new "A#_rep#" naming convention
        dataFilePattern = fullfile(runPath, '*.csv');
        dataFiles = dir(dataFilePattern);

        % Filter out files that start with '.'
        dataFiles = dataFiles(arrayfun(@(x) x.name(1) ~= '.', dataFiles));
        
        if isempty(dataFiles)
            warning('No data files found in: %s', runPath);
            continue;
        end
        
        for j = 1:length(dataFiles)
            filePath = fullfile(runPath, dataFiles(j).name);
            fileName = dataFiles(j).name;
            
            % Extract the cell name (A#) from the filename using a regular expression
            nameMatch = regexp(fileName, '([A-Za-z]\d+)_rep\d+\.csv$', 'tokens');
            
            % If the filename matches the pattern, extract the cell name and add to the data table
            if ~isempty(nameMatch)
                cellName = nameMatch{1}{1};
                
                % Add the data to the rows
                dataRows(end+1, :) = {cellName, runName, concentration, filePath}; %#ok<AGROW>
            end
        end
    end
end

% Create a table from the collected data rows
dataTable = cell2table(dataRows, ...
    'VariableNames', {'CellName', 'RunName', 'Concentration', 'FilePath'});
end
% Author: Lina M. Acosta Perez
%Description: Creates a single plot for each measurement run, showing the average current vs. voltage for all concentrations within that run. It also saves the summary data to an Excel file.
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
    legendHandles = gobjects(0); % collect plotted line handles
    legendEntries = {};

    concentrations = unique(runData.Concentration);
    cmap = hsv(length(concentrations));

    summaryTableData = {'VGS'};
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

            v = data(:, 3);  % Voltage
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

        % Plot with standard deviation shading (exclude fill from legend)
        hFill = fill([voltages, fliplr(voltages)], ...
                [meanI + stdI, fliplr(meanI - stdI)], ...
                cmap(c, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        set(hFill, 'HandleVisibility', 'off');
        hLine = plot(voltages, meanI, 'Color', cmap(c, :), 'LineWidth', 2);
        legendHandles(end+1) = hLine;
        legendEntries{end+1} = sprintf('[%.2f]', conc);

        % Append to summary table
        if c == 1
            summaryTableData(2:length(voltages)+1, 1) = num2cell(voltages');
        end
        summaryTableData{1, end+1} = sprintf('[%.2f] Avg Current [A]', conc);
        summaryTableData{1, end+1} = sprintf('[%.2f] Std Dev [A]', conc);
        summaryTableData(2:length(meanI)+1, end-1) = num2cell(meanI');
        summaryTableData(2:length(stdI)+1, end) = num2cell(stdI');
    end

    xlabel('Gate Voltage [V]');
    ylabel('Average Current [A]');
    title(sprintf('Average Current by Concentration - Run: %s', runName));
    legend(legendHandles, legendEntries, 'Location', 'eastoutside');
    set(gca, 'Position', [0.1, 0.1, 0.65, 0.8]);
    grid on;
    hold off;
    ylim([minAbsoluteCurrent*0.90 maxAbsoluteCurrent*1.10]); % Set y-axis limit based on max absolute current

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

% Author: Lina M. Acosta Perez
% Description: Averages the data from all runs for a given cell, then plots the average current vs. concentration with error bars.
% Date: Sept. 12 2025

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
            if size(data, 2) < 4
                continue;
            end

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

function plot_concentration_vs_current_by_group_letter_logscale(dataTable, analysisFolder, targetVgsStr, max_IDS_raw, min_IDS_raw)

    if ~exist(analysisFolder, 'dir')
        mkdir(analysisFolder);
    end

    targetVGS = str2double(targetVgsStr);

    if isnan(targetVGS)
        error('Invalid VGS value entered.');
        return;
    end
    
    % Prepare Group Letters
    dataTable.GroupLetter = cellfun(@(s) regexp(s, '^[A-Za-z]', 'match', 'once'), ...
                                dataTable.CellName, 'UniformOutput', false);
    groupLetters = unique(dataTable.GroupLetter);

    % Loop through each Group
    for g = 1:length(groupLetters)
        group = groupLetters{g};
        groupData = dataTable(strcmp(dataTable.GroupLetter, group), :);

        concentrations = unique(groupData.Concentration);
        
        % Initialize arrays to hold summary data for the final plot
        summaryConcs = [];
        summaryCurrents = [];
        summaryStds = [];
        
        % Loop through each Concentration
        for c = 1:length(concentrations)
            conc = concentrations(c);
            files = groupData(groupData.Concentration == conc, :).FilePath;

            voltages = [];
            currents = []; % All sweeps for this concentration

            % Read all sweeps for this concentration
            for j = 1:height(files)
                data = readmatrix(files{j});
                if size(data, 2) < 4
                    continue;
                end

                v = data(:, 3); % Gate Voltage
                i = data(:, 4); % Current

                if isempty(voltages)
                    voltages = v(:)';
                end

                if length(v) == length(voltages)
                    currents(end+1, :) = i(:)';
                else
                    warning('Skipping mismatched file: %s', files{j});
                end
            end

            if isempty(currents)
                continue;
            end
            
            % Find the index closest to the target VGS
            [~, idx] = min(abs(voltages - targetVGS));
            
            % Extract the current values at the target VGS for all sweeps
            currents_at_targetVGS = currents(:, idx);

            % Calculate Mean and Standard Deviation of current at that single point
            meanI_at_targetVGS = mean(currents_at_targetVGS, 'omitnan');
            stdI_at_targetVGS = std(currents_at_targetVGS, 0, 'omitnan');
            
            % Store the concentration and the extracted current/std
            summaryConcs(end+1) = conc;
            summaryCurrents(end+1) = meanI_at_targetVGS;
            summaryStds(end+1) = stdI_at_targetVGS;
        end % End concentration loop

        % --- Final Plot Generation ---
        
        % Sort data for plotting (optional, but good practice for log plots)
        [summaryConcs, sortIdx] = sort(summaryConcs);
        summaryCurrents = summaryCurrents(sortIdx);
        summaryStds = summaryStds(sortIdx);
        
        figure;
        % Use errorbar to plot Concentration (X) vs. Current (Y)
        errorbar(summaryConcs, summaryCurrents, summaryStds, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6, 'CapSize', 8);

        % Set the X-axis to logarithmic scale (CRITICAL change for this function name)
        set(gca, 'XScale', 'log'); 
        
        % Set Correct Axis Labels
        xlabel('Concentration [pM] (Log Scale)');
        ylabel(sprintf('Avg Current at V_{GS} = %.2f V [A]', targetVGS));
        title(sprintf('Concentration vs Current - Group: %s', group));
        
        ylim([min_IDS_raw*0.90 max_IDS_raw*1.10]); 
        grid on;

        % Save the figure
        plotFile = fullfile(analysisFolder, sprintf('Group_%s_conc_vs_current_logscale.png', group));
        saveas(gcf, plotFile);
        fprintf('✅ Saved Concentration vs Current (Log Scale) plot for group %s: %s\n', group, plotFile);

        % Save the Excel file (Concentration, Mean Current, Std Dev)
        excelFile = fullfile(analysisFolder, sprintf('Group_%s_conc_vs_current_logscale.xlsx', group));
        summaryTable = table(summaryConcs(:), summaryCurrents(:), summaryStds(:), ...
            'VariableNames', {'Concentration', 'AvgCurrent', 'StdDev'});
        writetable(summaryTable, excelFile);
        fprintf('✅ Saved Excel for group %s: %s\n', group, excelFile);
    end % End group loop

end

% Author: Lina M. Acosta Perez
% Description: Plots average VGS vs. concentration for each cell and run at a specific targetCurrentStr (e.g., '1e-6'). 
% Date: Sept. 12 2025

function plot_concentration_vs_VTH(dataTable, analysisFolder, targetCurrentStr, max_VGS_raw, min_VGS_raw)
% plot_concentration_vs_VTH Generates plots and saves data for concentration vs. VGS at a specific current.
% This function uses a target current value provided as an input argument.
% For each concentration within each cell and run, it finds the average VGS
% at which the current is closest to the target current. It then plots
% Concentration vs. Average VGS with error bars and saves the plot as a PNG
% and the summary data as an Excel file.
%
% Inputs:
%   dataTable        - A MATLAB table containing the gathered data,
%                      expected to have columns like 'CellName', 'RunName',
%                      'Concentration', and 'FilePath'.
%   analysisFolder   - The path to the folder where the generated plots
%                      and Excel files will be saved.
%   targetCurrentStr - A string representing the target current value (e.g., '1e-6').

    % Create analysis folder if it doesn't exist
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
                        warning('Skipping file %s: Not enough columns (expected at least 4 for VGS and Current).', files{j});
                        continue;
                    end
                    
                    v = data(:, 3); % VGS (Gate Voltage)
                    i = data(:, 4); % Current (Drain Current)

                    % Ensure current data is numeric and not empty
                    if isempty(i) || ~isnumeric(i)
                        warning('Skipping file %s: Current data is empty or not numeric.', files{j});
                        continue;
                    end
                    
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
            ylabel(sprintf('Avg V_{GS} at Current = %.2e A', targetCurrent));
            title(sprintf('Cell: %s | %s', cellName, runName));
            ylim([min_VGS_raw*0.90 max_VGS_raw*1.10]); % Set y-axis limit based on max VGS from all data

            grid on; % Add a grid for better readability
            
            % Construct file names including the Current value to avoid overwrites
            current_str_raw = sprintf('Current_%.3e', targetCurrent);
            % Replace characters that might be problematic in filenames (e.g., '.', '+', '-')
            current_str = strrep(current_str_raw, '.', 'p');
            current_str = strrep(current_str, '+', '');
            current_str = strrep(current_str, '-', 'm');
            
            summaryPlotFile = fullfile(analysisFolder, sprintf('%s_%s_concentration_vs_VGS_%s.png', cellName, runName, current_str));
            summaryExcelFile = fullfile(analysisFolder, sprintf('%s_%s_concentration_vs_VGS_%s.xlsx', cellName, runName, current_str));
            
            % Save the plot as a PNG image
            saveas(gcf, summaryPlotFile);
            
            % Close the figure to prevent too many open figures
            close(gcf);
            
            % Create a table from summary data and save it to an Excel file
            summaryTable = table(summaryConcs(:), summaryVGSValues(:), summaryVGSStds(:), ...
                                 'VariableNames', {'Concentration', 'AvgVGS', 'StdDevVGS'});
            writetable(summaryTable, summaryExcelFile);
        end
    end
    fprintf('Finished plotting for the specified target current.\n');
end

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
        legendHandles = gobjects(0);
        legendEntries = {};
        avgData = [];

        for c = 1:length(concentrations)
            conc = concentrations(c);
            files = runData(runData.Concentration == conc, :).FilePath;

            voltages = [];
            currents = [];

            for j = 1:height(files)
                data = readmatrix(files{j});
                if size(data, 2) < 4
                    continue;
                end

                v = data(:, 3);
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

            if isempty(currents)
                continue;
            end

            meanI = mean(currents, 1, 'omitnan');
            stdI = std(currents, 0, 1, 'omitnan');

            % Plot with standard deviation shading and exclude fill from legend
            hFill = fill([voltages, fliplr(voltages)], ...
                    [meanI + stdI, fliplr(meanI - stdI)], ...
                    cmap(c, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            set(hFill, 'HandleVisibility', 'off');
            hLine = plot(voltages, meanI, 'Color', cmap(c, :), 'LineWidth', 2);
            legendHandles(end+1) = hLine;
            legendEntries{end+1} = sprintf('[%.2f]', conc);

            % Append to summary table
            if c == 1
                avgData(1).voltages = voltages';
                avgData(1).concentration = conc;
            end
            avgData(c).allCurrents = currents;
            avgData(c).concentration = conc;
            avgData(c).voltage = voltages;
            avgData(c).meanCurrent = meanI;
            avgData(c).stdCurrent = stdI;
        end

        xlabel('Gate Voltage [V]');
        ylabel('Average Current [A]');
        title(sprintf('Cell: %s | %s', cellName, runName));
        legend(legendHandles, legendEntries, 'Location', 'eastoutside');
        set(gca, 'Position', [0.1, 0.1, 0.65, 0.8]);
        grid on;
        hold off;
        ylim([minAbsoluteCurrent*0.90 maxAbsoluteCurrent*1.10]);

        plotFile = fullfile(analysisFolder, sprintf('%s_%s_avg_plot.png', cellName, runName));
        saveas(gcf, plotFile);

        excelFile = fullfile(analysisFolder, sprintf('%s_%s_avg_data.xlsx', cellName, runName));
        dataMatrix = avgData(1).voltage(:);
        headers = {'VGS'};

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
                if size(data, 2) < 4
                    continue;
                end

                % Modified: Change column 3 to 2 for x-axis data
                v = data(:, 3);  % Gate Voltage
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
        xlabel('Gate Voltage (V)');
        ylabel('Average Current [A]');
        title(sprintf('Group %s - Sweep Plot', group));
        legend(h, legendEntries, 'Location', 'eastoutside');
        set(gca, 'Position', [0.1, 0.1, 0.65, 0.8]);
        grid on;
        hold off;
        ylim([minAbsoluteCurrent*0.90 maxAbsoluteCurrent*1.10]);

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


    % --- Helper: natural / numeric-aware sort for cellstr run names
    function sortedNames = sort_nat(names)
        % Accept cell array of char or string array
        if isstring(names)
            names = cellstr(names);
        elseif ischar(names)
            names = cellstr(names);
        end

        n = numel(names);
        if n == 0
            sortedNames = names;
            return;
        end

        % Try full numeric conversion (all entries numeric)
        nums = str2double(names);
        if all(~isnan(nums))
            [~, idx] = sort(nums);
            sortedNames = names(idx);
            return;
        end

        % Extract the first numeric token from each string (if any)
        numToken = nan(n,1);
        for k = 1:n
            t = regexp(names{k}, '(\d+)', 'tokens', 'once');
            if ~isempty(t)
                numToken(k) = str2double(t{1});
            end
        end

        % Use numeric token as primary key (missing -> Inf), lexicographic as secondary
        primaryKey = numToken;
        primaryKey(isnan(primaryKey)) = Inf;

        % secondary key: lexicographic
        [~, lexOrder] = sort(names);

        % sort by primaryKey, then by lexOrder
        [~, order] = sortrows([primaryKey, lexOrder]);
        sortedNames = names(order);
    end
end 

% Author: Lina M. Acosta Perez
% Description: Plots average current vs. concentration for each cell and run. The current is taken at a targetVgsStr (e.g., '0.5'). 
% Date: Sept. 12 2025

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
                    v = data(:, 3);
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
            ylim([min_IDS_raw*0.90 max_IDS_raw*1.10]);
         
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

                    v = data(:, 3);
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
            ylim([min_IDS_raw*0.90 max_IDS_raw*1.10]);
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

                v = data(:, 3);
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
        ylim([min_IDS_raw*0.90 max_IDS_raw*1.10]);
      
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
% Description: Plots average current versus run number for each concentration at a specific target voltage (VGS). It generates both combined and individual plots, as well as Excel tables.
% Date: Sept. 12 2025

function plot_run_vs_current_by_concentration(dataTable, analysisFolder, targetVgsStr, max_IDS_raw, min_IDS_raw)
    
    if ~exist(analysisFolder, 'dir')
        mkdir(analysisFolder);
    end

    targetVGS = str2double(targetVgsStr);

    if isnan(targetVGS)
        error('Invalid VGS value entered.');
        return;
    end

    runNames = unique(dataTable.RunName);
    runNames = sort(runNames); % Sort run names numerically if applicable
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
        ylim([min_IDS_raw*0.90 max_IDS_raw*1.10]);
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

% Author: Lina M. Acosta Perez
% Description: Finds the average VGS at a specific targetCurrentStr for each concentration within a group (determined by the first letter of the CellName). It then plots concentration vs. average VGS with error bars. The x-axis is in logarithmic scale.
% Date: Sept. 12 2025

function plot_concentration_vs_VTH_by_letter_logscale(dataTable, analysisFolder, targetCurrentStr, max_VGS_raw, min_VGS_raw)
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
        ylim([min_VGS_raw*0.90 max_VGS_raw*1.10]);

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



% Author: Lina M. Acosta Perez
% Description: Finds the average VGS at a specific targetCurrentStr for each concentration within a group (determined by the first letter of the CellName). It then plots concentration vs. average VGS with error bars. 
% Date: Sept. 12 2025

function plot_concentration_vs_VTH_by_letter(dataTable, analysisFolder, targetCurrentStr, max_VGS_raw, min_VGS_raw)
% plot_concentration_vs_current_by_group_letter Generates plots and saves data for concentration vs. VGS at a specific current, grouped by letter.
% This function prompts the user for a target current value. For each
% concentration within each group letter, it finds the average VGS at which
% the current is closest to the target current. It then plots Concentration
% vs. Average VGS with error bars and saves the plot as a PNG and the
% summary data as an Excel file.
%
% Inputs:
%   dataTable        - A MATLAB table containing the gathered data,
%                      expected to have columns like 'CellName', 'Concentration',
%                      and 'FilePath'. This function adds a 'GroupLetter' column.
%   analysisFolder   - The path to the folder where the generated plots
%                      and Excel files will be saved.

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
        % The original code did not have 'xscale log;' here, so it's not added.
        ylabel(sprintf('Avg V_{GS} at Current = %.2e A', targetCurrent));
        title(sprintf('Group %s - Concentration vs V_{GS}', group));
        ylim([min_VGS_raw*0.90 max_VGS_raw*1.10]);

        grid on; % Add a grid for better readability

        % Construct file names including the Current value to avoid overwrites
        current_str_raw = sprintf('Current_%.3e', targetCurrent);

        summaryPlotFile = fullfile(analysisFolder, sprintf('Group_%s_concentration_vs_VGS_%s.png', group, current_str_raw));
        summaryExcelFile = fullfile(analysisFolder, sprintf('Group_%s_concentration_vs_VGS_%s.csv', group, current_str_raw)); % Original was .csv

        % Save the plot as a PNG image
        saveas(gcf, summaryPlotFile);
        % Close the figure to prevent too many open figures
        close(gcf);

        % Create a table from summary data and save it to an Excel file
        summaryTable = table(summaryConcs(:), summaryVGSValues(:), summaryVGSStds(:), ...
                             'VariableNames', {'Concentration', 'AvgVGS', 'StdDevVGS'});
        writetable(summaryTable, summaryExcelFile);
    end
    fprintf('Finished plotting for the specified target current, grouped by letter.\n');
end

% Author: Lina M. Acosta Perez
% Description: Plots average VGS vs. concentration for each cell and run at a specific targetCurrentStr (e.g., '1e-6'). The x-axis (concentration) in a logarithmic scale.
% Date: Sept. 12 2025

function plot_concentration_vs_VTH_logscale(dataTable, analysisFolder, targetCurrentStr, max_VGS_raw, min_VGS_raw)
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
            ylim([min_VGS_raw*0.90 max_VGS_raw*1.10]);

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


% Author: Lina M. Acosta Perez
% Description: Plots average VGS vs. Run Name for each concentration at a specific targetCurrentStr.
% Date: Sept. 12 2025

function plot_run_vs_VTH_by_concentration(dataTable, analysisFolder, targetCurrentStr, max_VGS_raw, min_VGS_raw)
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
        ylim([min_VGS_raw*0.90 max_VGS_raw*1.10]);
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

