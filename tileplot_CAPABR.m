function [fig, t] = tileplot_CAPABR(y_avg_sel, all_freqs, A_descending)
% Plot average CAP or ABR data in tiled plot format. 
%
% Input:
%   y_avg_sel - cell array, each cell entry is sorted stack of average traces for a frequency
%
% 8/18/2023 George Liu
% Dependencies: none

SAMPLE_PERIOD = 5.12; % sample period of 5.12 us/sample in single trace
SAMPLE_PERIOD_MS = SAMPLE_PERIOD/1000; % sample period in ms

% Obtain plotting parameters
n_freq = length(all_freqs);
n_levels = length(A_descending);
num_samples = size(y_avg_sel{1}, 2);
X = SAMPLE_PERIOD_MS * (1:num_samples); % time in ms

% Make sure bounds of plot are within y limits
ylim_max = [];
ymax_data = max(cell2mat(y_avg_sel), [], "all");
ymin_data = min(cell2mat(y_avg_sel), [], "all");
ylim_max(2) = ymax_data;
ylim_max(1) = ymin_data;

% Plot
fig = figure;
t = tiledlayout(n_levels, n_freq, 'TileIndexing', 'columnmajor');

for j = 1:n_freq
    for k = 1:n_levels
        nexttile(t)
        plot(X, y_avg_sel{j}(k, :), 'LineWidth', 3, 'Color', [0 0.4470 0.7410])
        ylim(ylim_max)

        % show ylabels for first column only
        if j == 1
            ylabel([num2str(A_descending(k)), ' dB (nV)'], 'FontSize', 24)
            % rotate y label to make horizontal
            hYLabel = get(gca,'YLabel');
            set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right', 'FontSize', 24)
        end

        % Adjust plot appearance
        % Remove extraneous axis tick marks and x-axis from all but bottom
        % tile
        set(gca,'box','off')
        if mod(k, n_levels) ~= 0 % 
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'XColor','none')
        end

        % Remove y-axis labels from all but first column
        if j > 1
            set(gca,'yticklabel', [])
        end

        % Set text size of current axis
        set(gca,'FontSize',24)

        % Show frequency as title above top-most tile of each column
        if k==1
            title([num2str(all_freqs(j)), ' Hz']) 
        end
    end
end

% Format figure
xlabel(t, 'Time (ms)', 'FontSize', 24)
%     ylabel(t, 'Voltage (nV)', 'FontSize', 24)

end