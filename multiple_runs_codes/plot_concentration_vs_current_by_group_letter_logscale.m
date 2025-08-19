function plot_concentration_vs_current_by_group_letter_logscale(dataTable, analysisFolder, targetVgsStr)
   
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
        xscale log;
        ylabel(sprintf('Avg Current at V_{GS} = %.2f V', targetVGS));
        ylim([0 0.0035]);
       
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
