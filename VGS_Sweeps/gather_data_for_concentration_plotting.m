function dataTable = gather_data_for_concentration_plotting(mainFolderPath)
    % Initializes storage for data as a struct array for robustness.
    tempData = struct('CellName', {}, 'RunName', {}, 'Parameters', {}, 'Concentration', {}, 'FilePath', {});

    % Get concentration folders
    concentrationFolders = dir(mainFolderPath);
    concentrationFolders = concentrationFolders([concentrationFolders.isdir]);
    concentrationFolders = concentrationFolders(~ismember({concentrationFolders.name}, {'.', '..'}));
    concentrationFolders = concentrationFolders(arrayfun(@(x) x.name(1) ~= '.', concentrationFolders));

    if isempty(concentrationFolders)
        warning('No concentration folders found in the specified path: %s', mainFolderPath);
        dataTable = table();
        return;
    end
    
    for c = 1:length(concentrationFolders)
        concName = concentrationFolders(c).name;
        concPath = fullfile(mainFolderPath, concName);

        % Extract numeric concentration from folder name [X]
        concMatch = regexp(concName, '\[(\d+\.?\d*)\]', 'tokens');
        if isempty(concMatch)
            warning('Could not extract concentration from folder name: %s', concName);
            continue;
        end
        concentration = str2double(concMatch{1}{1});

        % Get run folders inside the concentration folder
        runFolders = dir(concPath);
        runFolders = runFolders([runFolders.isdir]);
        runFolders = runFolders(~ismember({runFolders.name}, {'.', '..'}));
        runFolders = runFolders(arrayfun(@(x) x.name(1) ~= '.', runFolders));
        
        if isempty(runFolders)
            warning('No run folders found in concentration folder: %s', concPath);
            continue;
        end

        for r = 1:length(runFolders)
            runName = runFolders(r).name;
            runPath = fullfile(concPath, runName);

            % Get parameter folders inside the run folder
            parameterFolders = dir(runPath);
            parameterFolders = parameterFolders([parameterFolders.isdir]);
            parameterFolders = parameterFolders(~ismember({parameterFolders.name}, {'.', '..'}));
            parameterFolders = parameterFolders(arrayfun(@(x) x.name(1) ~= '.', parameterFolders));

            if isempty(parameterFolders)
                warning('No parameter folders found in run folder: %s', runPath);
                continue;
            end
            
            for p = 1:length(parameterFolders)
                paramsName = parameterFolders(p).name;
                paramsPath = fullfile(runPath, paramsName);

                % Get cell folders inside the parameter folder
                cellFolders = dir(paramsPath);
                cellFolders = cellFolders([cellFolders.isdir]);
                cellFolders = cellFolders(~ismember({cellFolders.name}, {'.', '..'}));
                cellFolders = cellFolders(arrayfun(@(x) x.name(1) ~= '.', cellFolders));
                
                if isempty(cellFolders)
                    warning('No cell folders found in parameter folder: %s', paramsPath);
                    continue;
                end

                for cellIdx = 1:length(cellFolders)
                    cellName = cellFolders(cellIdx).name;
                    cellPath = fullfile(paramsPath, cellName);

                    % Find all CSV files inside the cell folder
                    csvFiles = dir(fullfile(cellPath, '*.csv'));
                    csvFiles = csvFiles(arrayfun(@(x) x.name(1) ~= '.', csvFiles));

                    if isempty(csvFiles)
                        warning('No CSV files found in: %s', cellPath);
                        continue;
                    end

                    for fileIdx = 1:length(csvFiles)
                        filePath = fullfile(cellPath, csvFiles(fileIdx).name);
                        
                        % Add a new struct to our temporary array
                        newRow.CellName = string(cellName);
                        newRow.RunName = string(runName);
                        newRow.Parameters = string(paramsName);
                        newRow.Concentration = concentration;
                        newRow.FilePath = string(filePath);
                        
                        tempData(end+1) = newRow;
                    end
                end
            end
        end
    end
    
    % Create the final table from the struct array
    if isempty(tempData)
        dataTable = table('Size', [0, 5], 'VariableTypes', {'string', 'string', 'string', 'double', 'string'}, ...
            'VariableNames', {'CellName', 'RunName', 'Parameters', 'Concentration', 'FilePath'});
        warning('No data found in the specified folder structure.');
    else
        dataTable = struct2table(tempData);
    end
end
