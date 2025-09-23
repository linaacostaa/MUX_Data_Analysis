% Author: Lina M. Acosta Perez
% Description: This version a current you enter, finds the closest data points, and plots the average voltage against concentration.
% The plots and data are also saved. It now correctly handles standard deviation data.
% Date: Sept. 13 2025
function plot_avg_from_multiple_csvs_current_sampling_2()
    % This script processes data from multiple CSV files to analyze and visualize
    % semiconductor device characteristics. It is adapted to handle CSVs where
    % data is organized in pairs of columns: Average Current and Standard Deviation,
    % for each concentration. It prompts for a current value and plots corresponding voltages.
    
    % Ask the user to select multiple CSV files
    [files, path] = uigetfile('*.csv', 'Select CSV Files', 'MultiSelect', 'on');
    if isequal(files, 0)
        error('No files selected.');
    end
    if ischar(files)
        files = {files}; % convert to cell array if only one file
    end
    
    allData = containers.Map('KeyType', 'char', 'ValueType', 'any');
    allStdData = containers.Map('KeyType', 'char', 'ValueType', 'any');
    vgs_values = [];
    
    % Process each CSV file
    for f = 1:length(files)
        filePath = fullfile(path, files{f});
        T = readcell(filePath);
        
        % Extract VGS column, assuming it's the first
        vgs_values = cell2mat(T(2:end, 1));
        
        % Loop through headers in pairs (Avg and Std)
        for col = 2:2:size(T, 2)
            avgHeader = T{1, col};
            stdHeader = T{1, col+1};
            
            match = regexp(avgHeader, '\[(.*?)\]', 'tokens');
            if isempty(match), continue; end
            concStr = match{1}{1};
            
            current_data = cell2mat(T(2:end, col));
            std_data = cell2mat(T(2:end, col+1));
            
            % For combined data
            if ~isKey(allData, concStr)
                allData(concStr) = {};
                allStdData(concStr) = {};
            end
            allData(concStr) = [allData(concStr), {current_data}];
            allStdData(concStr) = [allStdData(concStr), {std_data}];
        end
    end
    
    %% Ask for current to extract voltage at
    prompt = {'Enter the Current value to extract voltage at:'};
    dlgtitle = 'Current Selection';
    dims = [1 50];
    definput = {'1e-6'};
    answer = inputdlg(prompt, dlgtitle, dims, definput);
    if isempty(answer)
        return;
    end
    
    targetCurrent = str2double(answer{1});
    if isnan(targetCurrent)
        errordlg('Invalid Current value entered.');
        return;
    end
    
    concKeys = sort(keys(allData));
    concVals = str2double(concKeys);
    
    voltagesAtCurrent = zeros(size(concVals));
    stdsOfVoltage = zeros(size(concVals));
    
    for i = 1:length(concKeys)
        conc = concKeys{i};
        currents_per_file = allData(conc);
        
        % Initialize lists for voltage values at the target current
        voltageList = [];
        
        for k = 1:length(currents_per_file)
            currents = currents_per_file{k};
            [~, idx] = min(abs(currents - targetCurrent));
            voltageList = [voltageList, vgs_values(idx)];
        end
        
        voltagesAtCurrent(i) = mean(voltageList, 'omitnan');
        stdsOfVoltage(i) = std(voltageList, 'omitnan');
    end
    
    [sortedConc, sortIdx] = sort(concVals);
    sortedVoltages = voltagesAtCurrent(sortIdx);
    sortedStds = stdsOfVoltage(sortIdx);
    
    %% Plot combined data
    figure;
    errorbar(sortedConc, sortedVoltages, sortedStds, '-o', 'LineWidth', 2);
    xlabel('Concentration');
    ylabel(sprintf('Average Voltage at I = %.2e A', targetCurrent));
    title(sprintf('Concentration vs Voltage at I = %.2e A', targetCurrent));
    grid on;
    
    % Set y-axis limits BEFORE saving the figure
    ylim([0.5 1.2]);
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    saveas(gcf, fullfile(path, ['concentration_vs_voltage_' timestamp '.png']));
    
    figure;
    errorbar(sortedConc, sortedVoltages, sortedStds, '-o', 'LineWidth', 2);
    set(gca, 'XScale', 'log');
    xlabel('Concentration (log scale)');
    ylabel(sprintf('Average Voltage at I = %.2e A', targetCurrent));
    title(sprintf('Concentration vs Voltage (log) at I = %.2e A', targetCurrent));
    grid on;
    
    % Set y-axis limits BEFORE saving the figure
    ylim([0.5 1.2]);
    saveas(gcf, fullfile(path, ['concentration_vs_voltage_log_' timestamp '.png']));
    
    summaryTable = table(sortedConc(:), sortedVoltages(:), sortedStds(:), ...
        'VariableNames', {'Concentration', 'AvgVoltage', 'StdDev'});
    writetable(summaryTable, fullfile(path, ['concentration_vs_voltage_' timestamp '.xlsx']));
    
    %% Plot each experiment separately at a given current
    numExperiments = length(files);
    
    figure;
    hold on;
    cmap_exp = hsv(numExperiments);
    
    for i = 1:numExperiments
        filePath = fullfile(path, files{i});
        T_exp = readcell(filePath);
        vgs_exp = cell2mat(T_exp(2:end, 1));
        
        voltagesAtCurrentExp = NaN(size(concVals));
        stdsAtCurrentExp = NaN(size(concVals));
        
        for col = 2:2:size(T_exp, 2)
            avgHeader = T_exp{1, col};
            match = regexp(avgHeader, '\[(.*?)\]', 'tokens');
            if isempty(match), continue; end
            
            concStr_exp = match{1}{1};
            
            currents = cell2mat(T_exp(2:end, col));
            stds = cell2mat(T_exp(2:end, col+1));
            
            [~, idx] = min(abs(currents - targetCurrent));
            
            fullConcIdx = find(concVals == str2double(concStr_exp));
            
            if ~isempty(fullConcIdx)
                voltagesAtCurrentExp(fullConcIdx) = vgs_exp(idx);
                stdsAtCurrentExp(fullConcIdx) = stds(idx);
            end
        end
        
        errorbar(concVals, voltagesAtCurrentExp, stdsAtCurrentExp, '-o', 'LineWidth', 2, 'Color', cmap_exp(i, :));
    end
    
    hold off;
    legend(files, 'Location', 'eastoutside');
    xlabel('Concentration');
    ylabel(sprintf('Voltage at I = %.2e A', targetCurrent));
    title('Concentration vs Voltage by Experiment');
    grid on;
    
    % Set y-axis limits BEFORE saving the figure
    ylim([0.5 1.2]);
    saveas(gcf, fullfile(path, ['experiment_overlap_voltage_plot_' timestamp '.png']));
    
    figure;
    hold on;
    for i = 1:numExperiments
        filePath = fullfile(path, files{i});
        T_exp = readcell(filePath);
        vgs_exp = cell2mat(T_exp(2:end, 1));
        
        voltagesAtCurrentExp = NaN(size(concVals));
        stdsAtCurrentExp = NaN(size(concVals));
        
        for col = 2:2:size(T_exp, 2)
            avgHeader = T_exp{1, col};
            match = regexp(avgHeader, '\[(.*?)\]', 'tokens');
            if isempty(match), continue; end
            
            concStr_exp = match{1}{1};
            currents = cell2mat(T_exp(2:end, col));
            stds = cell2mat(T_exp(2:end, col+1));
            
            [~, idx] = min(abs(currents - targetCurrent));
            
            fullConcIdx = find(concVals == str2double(concStr_exp));
            
            if ~isempty(fullConcIdx)
                voltagesAtCurrentExp(fullConcIdx) = vgs_exp(idx);
                stdsAtCurrentExp(fullConcIdx) = stds(idx);
            end
        end
        
        errorbar(concVals, voltagesAtCurrentExp, stdsAtCurrentExp, '-o', 'LineWidth', 2, 'Color', cmap_exp(i, :));
    end
    set(gca, 'XScale', 'log');
    hold off;
    legend(files, 'Location', 'eastoutside');
    xlabel('Concentration (log scale)');
    ylabel(sprintf('Voltage at I = %.2e A', targetCurrent));
    title('Concentration vs Voltage by Experiment (Log Scale)');
    grid on;
    
    % Set y-axis limits BEFORE saving the figure
    ylim([0.5 1.2]);
    saveas(gcf, fullfile(path, ['experiment_overlap_voltage_log_plot_' timestamp '.png']));
    
    %% Export the data for each experiment at the selected current value
    exportExpData = table(concVals(:), 'VariableNames', {'Concentration'});
    for i = 1:numExperiments
        filePath = fullfile(path, files{i});
        T_exp = readcell(filePath);
        vgs_exp = cell2mat(T_exp(2:end, 1));
        
        voltagesAtCurrentExp = NaN(size(concVals));
        stdsAtCurrentExp = NaN(size(concVals));
        
        for col = 2:2:size(T_exp, 2)
            avgHeader = T_exp{1, col};
            match = regexp(avgHeader, '\[(.*?)\]', 'tokens');
            if isempty(match), continue; end
            
            concStr_exp = match{1}{1};
            currents = cell2mat(T_exp(2:end, col));
            stds = cell2mat(T_exp(2:end, col+1));
            
            [~, idx] = min(abs(currents - targetCurrent));
            
            fullConcIdx = find(concVals == str2double(concStr_exp));
            
            if ~isempty(fullConcIdx)
                voltagesAtCurrentExp(fullConcIdx) = vgs_exp(idx);
                stdsAtCurrentExp(fullConcIdx) = stds(idx);
            end
        end
        
        avgVarName = sprintf('AvgVoltage_%s', strrep(files{i}, '.csv', ''));
        stdVarName = sprintf('StdDev_%s', strrep(files{i}, '.csv', ''));
        exportExpData.(avgVarName) = voltagesAtCurrentExp';
        exportExpData.(stdVarName) = stdsAtCurrentExp';
    end
    writetable(exportExpData, fullfile(path, ['experiment_data_at_current_' timestamp '.xlsx']));
end