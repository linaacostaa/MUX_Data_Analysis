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
concFolderPath = uigetdir('', 'Select the folder containing [Concentration]/Run folders');
if concFolderPath == 0
    error('No folder selected for concentration analysis. Exiting...');
end
%% Step 3: Gather data for concentration plots
dataTable = gather_data_for_concentration_plotting(concFolderPath);
%% Step 4: Create a folder to save concentration plots inside the selected folder
analysisFolder_concentration = fullfile(concFolderPath, '4_voltage_vs_current');
if ~exist(analysisFolder_concentration, 'dir')
    mkdir(analysisFolder_concentration);
end
%% Step 4.1: Generate average current plots by concentration - only averages the 5 sweeps of a cell name at one run
plot_voltage_vs_current(dataTable, analysisFolder_concentration);
%% Step 7: Create a folder to save concentration plots inside the selected folder
analysisFolder_volt_vs_avg_current_letter = fullfile(concFolderPath, '7_volt_vs_current_avg_letter');
if ~exist(analysisFolder_volt_vs_avg_current_letter, 'dir')
    mkdir(analysisFolder_volt_vs_avg_current_letter);
end
%% Step 7.1: Generate average current plots by concentration
plot_volt_vs_avg_current_by_group_letter(dataTable, analysisFolder_volt_vs_avg_current_letter);
%% Step 11: Create a folder to save concentration plots inside the selected folder
analysisFolder_avg_concentration = fullfile(concFolderPath, '11_sweep_avg_per_run');
if ~exist(analysisFolder_avg_concentration, 'dir')
    mkdir(analysisFolder_avg_concentration);
end
%% Step 11.1: Generate average of all the sweeps for a given run number
plot_avg_sweeps_by_run(dataTable, analysisFolder_avg_concentration);

%% Sampling Method Selection

samplingChoice = questdlg('Which sampling method would you like to perform?', ...
    'Sampling Method', ...
    'IDS sampling', 'VGS sampling', 'IDS sampling');

idsValue = '';
vgsValue = '';
%% IDS sampling

switch samplingChoice
    case 'IDS sampling'
        idsValue = inputdlg('Enter a value for IDS sampling:', 'IDS Value');
        if isempty(idsValue)
            error('IDS value not provided. Exiting...');
        end
        idsValue = idsValue{1}; % Extract value from cell array
        
        % IDS Sampling functions
        analysisFolder_concentration_vs_VTH = fullfile(concFolderPath, ['12_conc_vs_VTH_IDS_sampling_' idsValue]);
        if ~exist(analysisFolder_concentration_vs_VTH, 'dir')
            mkdir(analysisFolder_concentration_vs_VTH);
        end
        plot_concentration_vs_VTH(dataTable, analysisFolder_concentration_vs_VTH, idsValue);
        
        analysisFolder_concentration_vs_VTH_logscale = fullfile(concFolderPath, ['13_conc_vs_VTH_logscale_IDS_sampling_' idsValue]);
        if ~exist(analysisFolder_concentration_vs_VTH_logscale, 'dir')
            mkdir(analysisFolder_concentration_vs_VTH_logscale);
        end
        plot_concentration_vs_VTH_logscale(dataTable, analysisFolder_concentration_vs_VTH_logscale, idsValue);
        
        analysisFolder_conc_vs_avg_VTH_well = fullfile(concFolderPath, ['14_conc_vs_VTH_by_letter_IDS_sampling_' idsValue]);
        if ~exist(analysisFolder_conc_vs_avg_VTH_well, 'dir')
            mkdir(analysisFolder_conc_vs_avg_VTH_well);
        end
        plot_concentration_vs_VTH_by_letter(dataTable, analysisFolder_conc_vs_avg_VTH_well, idsValue);
        
        analysisFolder_conc_vs_avg_VTH_well_logscale = fullfile(concFolderPath, ['15_conc_vs_VTH_by_letter_log_IDS_sampling_' idsValue]);
        if ~exist(analysisFolder_conc_vs_avg_VTH_well_logscale, 'dir')
            mkdir(analysisFolder_conc_vs_avg_VTH_well_logscale);
        end
        plot_concentration_vs_VTH_by_letter_logscale(dataTable, analysisFolder_conc_vs_avg_VTH_well_logscale, idsValue);
        
        analysisFolder_plot_run_vs_VTH_by_concentration = fullfile(concFolderPath, ['16_run_vs_VTH_IDS_sampling_' idsValue]);
        if ~exist(analysisFolder_plot_run_vs_VTH_by_concentration, 'dir')
            mkdir(analysisFolder_plot_run_vs_VTH_by_concentration);
        end
        plot_run_vs_VTH_by_concentration(dataTable, analysisFolder_plot_run_vs_VTH_by_concentration, idsValue);
        
    %% VGS sampling
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
        plot_concentration_vs_current(dataTable, analysisFolder_concentration_vs_current, vgsValue);
         
        analysisFolder_concentration_vs_current_logscale = fullfile(concFolderPath, ['6_conc_vs_current_logscale_VGS_sampling_' vgsValue]);
        if ~exist(analysisFolder_concentration_vs_current_logscale, 'dir')
            mkdir(analysisFolder_concentration_vs_current_logscale);
        end
        plot_concentration_vs_current_logscale(dataTable, analysisFolder_concentration_vs_current_logscale, vgsValue);
       
        analysisFolder_conc_vs_avg_current_well = fullfile(concFolderPath, ['8_conc_vs_current_avg_letter_VGS_sampling_' vgsValue]);
        if ~exist(analysisFolder_conc_vs_avg_current_well, 'dir')
            mkdir(analysisFolder_conc_vs_avg_current_well);
        end
        plot_concentration_vs_current_by_group_letter(dataTable, analysisFolder_conc_vs_avg_current_well,vgsValue);
        
        analysisFolder_conc_vs_avg_current_well_logscale = fullfile(concFolderPath, ['9_conc_vs_current_avg_well_logscale_VGS_sampling_' vgsValue]);
        if ~exist(analysisFolder_conc_vs_avg_current_well_logscale, 'dir')
            mkdir(analysisFolder_conc_vs_avg_current_well_logscale);
        end
        plot_concentration_vs_current_by_group_letter_logscale(dataTable, analysisFolder_conc_vs_avg_current_well_logscale, vgsValue);
        
        analysisFolder_plot_run_vs_current_by_concentration = fullfile(concFolderPath, ['10_run_vs_IDS_VGS_sampling_' vgsValue]);
        if ~exist(analysisFolder_plot_run_vs_current_by_concentration, 'dir')
            mkdir(analysisFolder_plot_run_vs_current_by_concentration);
        end
        plot_run_vs_current_by_concentration(dataTable, analysisFolder_plot_run_vs_current_by_concentration, vgsValue);

         analysisFolder_concentration_vs_current_by_cell_name = fullfile(concFolderPath, ['20_conc_vs_current_per_cell_VGS_sampling_' vgsValue]);
        if ~exist(analysisFolder_concentration_vs_current_by_cell_name, 'dir')
            mkdir(analysisFolder_concentration_vs_current_by_cell_name);
        end
        plot_concentration_vs_current_by_cell_name(dataTable, analysisFolder_concentration_vs_current_by_cell_name, vgsValue);

        analysisFolder_perc_change_by_cell_name = fullfile(concFolderPath, ['21_conc_vs_current_per_change_VGS_sampling_' vgsValue]);
        if ~exist(analysisFolder_perc_change_by_cell_name, 'dir')
            mkdir(analysisFolder_perc_change_by_cell_name);
        end
        plot_perc_change_by_cell_name(dataTable, analysisFolder_perc_change_by_cell_name, vgsValue);
        
    otherwise
        disp('No sampling method selected or invalid choice. Exiting script.');
end
