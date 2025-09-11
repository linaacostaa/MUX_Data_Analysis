function dataTable = gather_data_for_concentration_plotting(mainFolderPath)
    dataRows = {};  % Initialize storage
    % Get concentration folders
    concentrationFolders = dir(mainFolderPath);
    concentrationFolders = concentrationFolders([concentrationFolders.isdir]);
    concentrationFolders = concentrationFolders(~ismember({concentrationFolders.name}, {'.', '..'}));
    
    % Filter out folders starting with '.'
    concentrationFolders = concentrationFolders(arrayfun(@(x) x.name(1) ~= '.', concentrationFolders));
    
    for c = 1:length(concentrationFolders)
        concName = concentrationFolders(c).name;
        concPath = fullfile(mainFolderPath, concName);
        
        % Extract numeric concentration value from folder name
        concentration = str2double(concName);
        if isnan(concentration), continue; end
        
        % Get group folders (e.g., VdsSweep_Vgs0.60V)
        groupFolders = dir(concPath);
        groupFolders = groupFolders([groupFolders.isdir]);
        groupFolders = groupFolders(~ismember({groupFolders.name}, {'.', '..'}));
        
        % Filter out folders starting with '.'
        groupFolders = groupFolders(arrayfun(@(x) x.name(1) ~= '.', groupFolders));
        
        for g = 1:length(groupFolders)
            groupName = groupFolders(g).name;
            groupPath = fullfile(concPath, groupName);
            
            % Get cell folders (e.g., A1, A2)
            cellFolders = dir(groupPath);
            cellFolders = cellFolders([cellFolders.isdir]);
            cellFolders = cellFolders(~ismember({cellFolders.name}, {'.', '..'}));
            
            % Filter out folders starting with '.'
            cellFolders = cellFolders(arrayfun(@(x) x.name(1) ~= '.', cellFolders));
            
            for i = 1:length(cellFolders)
                cellName = cellFolders(i).name;
                cellPath = fullfile(groupPath, cellName);
                
                % Find CSV files
                csvFiles = dir(fullfile(cellPath, '*.csv'));
                
                % Filter out CSV files starting with '.'
                csvFiles = csvFiles(arrayfun(@(x) x.name(1) ~= '.', csvFiles));
                
                for j = 1:length(csvFiles)
                    filePath = fullfile(cellPath, csvFiles(j).name);
                    dataRows(end+1, :) = {cellName, groupName, concentration, filePath}; %#ok<AGROW>
                end
            end
        end
    end
    
    dataTable = cell2table(dataRows, ...
        'VariableNames', {'CellName', 'GroupName', 'Concentration', 'FilePath'});
end