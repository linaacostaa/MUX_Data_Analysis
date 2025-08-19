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
        % Extract numeric concentration value from folder name [X]
        concMatch = regexp(concName, '\[(\d+\.?\d*)\]', 'tokens');
        if isempty(concMatch), continue; end
        concentration = str2double(concMatch{1}{1});
        % Get run folders inside this concentration folder
        runFolders = dir(concPath);
        runFolders = runFolders([runFolders.isdir]);
        runFolders = runFolders(~ismember({runFolders.name}, {'.', '..'}));
        
        % Filter out folders starting with '.'
        runFolders = runFolders(arrayfun(@(x) x.name(1) ~= '.', runFolders));
        
        for r = 1:length(runFolders)
            runName = runFolders(r).name;
            runPath = fullfile(concPath, runName);
            % Find JSON
            jsonFile = dir(fullfile(runPath, '*.json'));

            % Filter out JSON files starting with '.'
            jsonFile = jsonFile(arrayfun(@(x) x.name(1) ~= '.', jsonFile));

            if isempty(jsonFile)
                warning('No JSON file in: %s', runPath);
                continue;
            end
            try
                raw = fileread(fullfile(runPath, jsonFile.name));
                jsonData = jsondecode(raw);
                if ~isfield(jsonData, 'cell_names'), continue; end
            catch
                warning('Failed to read JSON in: %s', runPath);
                continue;
            end
            cellNames = jsonData.cell_names;
            for i = 1:length(cellNames)
                cellName = cellNames{i};
                % Find CSVs
                csvFiles = dir(fullfile(runPath, sprintf('*_%s_*.csv', cellName)));
                
                % Filter out CSV files that start with '.'
                csvFiles = csvFiles(arrayfun(@(x) x.name(1) ~= '.', csvFiles));
                
                for j = 1:length(csvFiles)
                    filePath = fullfile(runPath, csvFiles(j).name);
                    dataRows(end+1, :) = {cellName, runName, concentration, filePath}; %#ok<AGROW>
                end
            end
        end
    end
    dataTable = cell2table(dataRows, ...
        'VariableNames', {'CellName', 'RunName', 'Concentration', 'FilePath'});
end