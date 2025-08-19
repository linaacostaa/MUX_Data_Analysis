function plot_volt_vs_avg_current_by_group_letter(dataTable, analysisFolder)

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
        rawExport{1, 1} = 'VGS';
        colIdx = 2;

        for c = 1:length(concentrations)
            conc = concentrations(c);
            files = groupData(groupData.Concentration == conc, :).FilePath;

            voltages = [];
            currents = [];
            sweepLabels = {};

            for j = 1:height(files)
                data = readmatrix(files{j});
                if size(data, 2) < 4, continue; end

                v = data(:, 3);  % Voltage
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
        end

        xlabel('Gate Voltage [V]');
        ylabel('Average Current [A]');
        title(sprintf('Group %s - Sweep Plot', group));
        legend(h, legendEntries, 'Location', 'eastoutside');
        set(gca, 'Position', [0.1, 0.1, 0.65, 0.8]);
        grid on;
        hold off;

        % Save figure
        plotFile = fullfile(analysisFolder, sprintf('Group_%s_sweeps.png', group));
        saveas(gcf, plotFile);

        % Save raw data to CSV
        csvFile = fullfile(analysisFolder, sprintf('Group_%s_sweeps.csv', group));
        writecell(rawExport, csvFile);
    end
end
