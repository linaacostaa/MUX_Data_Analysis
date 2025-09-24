% Author: Lina M. Acosta Perez
% Description: The purpose of these scripts is to automate the visualization and summarization of electrical data from multiple measurement runs
% and cells. The scripts generate various plots and summary tables to help users quickly assess device performance across different concentrations and runs.
% Date: Sept. 12 2025

clc; clear; close all;
%% Step 1: Ask the user to select the scripts folder (where .m files are located)
scriptsFolder = uigetdir('', 'Select the Folder Containing .m Files');
if scriptsFolder == 0
    error('No scripts folder selected. Exiting...');
end
%% Step 2: Ask the user to select the folder with concentration folders
%concFolderPath = uigetdir('', 'Select the folder containing [Concentration]/Run folders');
%if concFolderPath == 0
 %   error('No folder selected for concentration analysis. Exiting...');
%end

parentFolderPath = uigetdir('', 'Select the parent folder that contains the [Concentration]/Run folders');
if parentFolderPath == 0
    error('No parent folder selected. Exiting...');
end

% Get all subfolders in the parent folder, which should be your concentration folders
subFolders = dir(parentFolderPath);
subFolders = subFolders([subFolders.isdir]);
subFolders = subFolders(~ismember({subFolders.name}, {'.', '..'}));

combinedDataRows = {};
for i = 1:length(subFolders)
    concFolderPath = fullfile(parentFolderPath, subFolders(i).name);
    fprintf('Processing folder: %s\n', concFolderPath);
    
    %% Gather data for concentration plots for the current folder
    dataTable = gather_data_for_concentration_plotting(concFolderPath);
    
    % Add a new column indicating the parent folder name for later comparisons
    if ~isempty(dataTable)
        dataTable.ParentFolder = repmat({subFolders(i).name}, height(dataTable), 1);
    end
    
    % Append this folder's data to the overall collection
    combinedDataRows = [combinedDataRows; table2cell(dataTable)];
end

% Recreate the combined data table.
if ~isempty(combinedDataRows)
    combinedDataTable = cell2table(combinedDataRows, ...
        'VariableNames', {'CellName', 'RunName', 'Concentration', 'FilePath', 'ParentFolder'});
else
    error('No data gathered from any subfolder.');
end

% Optional: Save the combined data table for future reference.
save(fullfile(parentFolderPath, 'combinedDataTable.mat'), 'combinedDataTable');

 %% Sampling Method Selection

    samplingChoice = questdlg('Which sampling method would you like to perform?', ...
        'Sampling Method', ...
        'IDS sampling', 'VGS sampling', 'IDS sampling');

    idsValue = '';
    vgsValue = '';% Example usage
    
% Assuming combinedDataTable is your table of numerical data
[I_V_min_raw, I_V_max_raw] = I_V_min_max(combinedDataTable);
fprintf('Current: min = %g, max = %g\n', I_V_min_raw, I_V_max_raw);

% For VDS sampling (when user has provided the target vgsValue or idsValue)
[targetValue] = vgsValue;  % or idsValue, depending on the sampling method
[VDS_min_raw, VDS_max_raw] = VDS_min_max(combinedDataTable, targetValue);
fprintf('Drain Voltage: min = %g, max = %g\n', VDS_min_raw, VDS_max_raw);

I_V_min_final = I_V_min_raw * 1.10; % Increase by 10%
I_V_max_final = I_V_max_raw * 1.10; % Increase by 10%
VDS_min_final = VDS_min_raw * 1.10; % Increase by 10%
VDS_max_final = VDS_min_raw * 1.10; % Increase by 10%

%I_V_min_final, I_V_max_final
%VDS_min_final, VDS_max_final



for i = 1:length(subFolders)
    %% Step 4.1: Generate average current plots by concentration - only averages the 5 sweeps of a cell name at one run
    plot_voltage_vs_current(dataTable, analysisFolder_concentration, I_V_min_final, I_V_max_final);
    %% Step 7: Create a folder to save concentration plots inside the selected folder
    analysisFolder_volt_vs_avg_current_letter = fullfile(concFolderPath, '7_volt_vs_current_avg_letter');
    if ~exist(analysisFolder_volt_vs_avg_current_letter, 'dir')
        mkdir(analysisFolder_volt_vs_avg_current_letter);
    end
    %% Step 7.1: Generate average current plots by concentration
    plot_volt_vs_avg_current_by_group_letter(dataTable, analysisFolder_volt_vs_avg_current_letter, I_V_min_final, I_V_max_final);
    %% Step 11: Create a folder to save concentration plots inside the selected folder
    analysisFolder_avg_concentration = fullfile(concFolderPath, '11_sweep_avg_per_run');
    if ~exist(analysisFolder_avg_concentration, 'dir')
        mkdir(analysisFolder_avg_concentration);
    end
    %% Step 11.1: Generate average of all the sweeps for a given run number
    plot_avg_sweeps_by_run(dataTable, analysisFolder_avg_concentration, I_V_min_final, I_V_max_final);

    %% Sampling Method Selection

    %samplingChoice = questdlg('Which sampling method would you like to perform?', ...
    %    'Sampling Method', ...
     %   'IDS sampling', 'VGS sampling', 'IDS sampling');

    %idsValue = '';
    %vgsValue = '';
    %% IDS sampling

    switch samplingChoice
        case 'VGS sampling'
            vgsValue = inputdlg('Enter a value for VGS sampling:', 'VGS Value');
            if isempty(vgsValue)
                error('VGS value not provided. Exiting...');
            end
            vgsValue = vgsValue{1}; % Extract value from cell array
            
            % VGS Sampling functions
            analysisFolder_concentration_vs_current = fullfile(concFolderPath, ['5_conc_vs_current_VGS_sampling_' vgsValue]);
            if ~exist(analysisFolder_concentration_vs_current, 'dir')
                mkdir(analysisFolder_concentration_vs_current);
            end
            plot_concentration_vs_current(dataTable, analysisFolder_concentration_vs_current, vgsValue, VDS_min_final, VDS_max_final);
            
            analysisFolder_concentration_vs_current_logscale = fullfile(concFolderPath, ['6_conc_vs_current_logscale_VGS_sampling_' vgsValue]);
            if ~exist(analysisFolder_concentration_vs_current_logscale, 'dir')
                mkdir(analysisFolder_concentration_vs_current_logscale);
            end
            plot_concentration_vs_current_logscale(dataTable, analysisFolder_concentration_vs_current_logscale, vgsValue, VDS_min_final, VDS_max_final);
        
            analysisFolder_conc_vs_avg_current_well = fullfile(concFolderPath, ['8_conc_vs_current_avg_letter_VGS_sampling_' vgsValue]);
            if ~exist(analysisFolder_conc_vs_avg_current_well, 'dir')
                mkdir(analysisFolder_conc_vs_avg_current_well);
            end
            plot_concentration_vs_current_by_group_letter(dataTable, analysisFolder_conc_vs_avg_current_well,vgsValue, VDS_min_final, VDS_max_final);
            
            analysisFolder_conc_vs_avg_current_well_logscale = fullfile(concFolderPath, ['9_conc_vs_current_avg_well_logscale_VGS_sampling_' vgsValue]);
            if ~exist(analysisFolder_conc_vs_avg_current_well_logscale, 'dir')
                mkdir(analysisFolder_conc_vs_avg_current_well_logscale);
            end
            plot_concentration_vs_current_by_group_letter_logscale(dataTable, analysisFolder_conc_vs_avg_current_well_logscale, vgsValue, VDS_min_final, VDS_max_final);
            
            analysisFolder_plot_run_vs_current_by_concentration = fullfile(concFolderPath, ['10_run_vs_IDS_VGS_sampling_' vgsValue]);
            if ~exist(analysisFolder_plot_run_vs_current_by_concentration, 'dir')
                mkdir(analysisFolder_plot_run_vs_current_by_concentration);
            end
            plot_run_vs_current_by_concentration(dataTable, analysisFolder_plot_run_vs_current_by_concentration, vgsValue, VDS_min_final, VDS_max_final);

            analysisFolder_concentration_vs_current_by_cell_name = fullfile(concFolderPath, ['20_conc_vs_current_per_cell_VGS_sampling_' vgsValue]);
            if ~exist(analysisFolder_concentration_vs_current_by_cell_name, 'dir')
                mkdir(analysisFolder_concentration_vs_current_by_cell_name);
            end
            plot_concentration_vs_current_by_cell_name(dataTable, analysisFolder_concentration_vs_current_by_cell_name, vgsValue, VDS_min_final, VDS_max_final);

            analysisFolder_perc_change_by_cell_name = fullfile(concFolderPath, ['21_conc_vs_current_per_change_V_sampling_' vgsValue]);
            if ~exist(analysisFolder_perc_change_by_cell_name, 'dir')
                mkdir(analysisFolder_perc_change_by_cell_name);
            end
            plot_perc_change_by_cell_name(dataTable, analysisFolder_perc_change_by_cell_name, vgsValue);
            
        otherwise
            disp('No sampling method selected or invalid choice. Exiting script.');
    end
end
