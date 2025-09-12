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
                nameMatch = regexp(fileName, '(A\d+)_rep\d+\.csv$', 'tokens');
                
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



