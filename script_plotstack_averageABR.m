% Script to plot average ABRs in CSV file in stacked plot
%
% 7/11/2021 George Liu
% Dependencies: import_averageABRcsv.m

%% Constants
SAMPLE_PERIOD = 40.96; % sample period of 40.96 us/sample in average trace

%% Load average trace data
[filename,path] = uigetfile('*.csv');
[M, A_csv, freq_csv] = import_averageABRcsv(filename, path);
num_samples = size(M, 2);
num_levels = size(M, 1);
plot_title = ['ABR @ ', num2str(freq_csv(1)), ' Hz'];

X = SAMPLE_PERIOD / 1000 * (1:num_samples); 
dB = cell(num_levels, 1);
dB(:) = {'dB (nV)'};
[A_descending, I] = sort(A_csv, 'descend');
M_sorted = M(I, :);
A_descending_cell = cellstr(num2str(A_descending));
newYlabels=cellfun(@(x,y) [x  ' ' y], A_descending_cell, dB, 'un', 0);

%% Plot ABRs in vertical stack
figure
s = stackedplot(X, M_sorted', 'Title', plot_title, 'DisplayLabels', newYlabels, 'LineWidth', 1.5);
% match y axis scale for all plots
ylim_max = [0, 0];
for j = 1:num_levels
    this_ylim_max = s.AxesProperties(j).YLimits;
    ylim_max(1) = min([ylim_max(1), this_ylim_max(1)]);
    ylim_max(2) = max([ylim_max(2), this_ylim_max(2)]);
end
% ylim_max = [-650, 650];
% ax = findobj(s.NodeChildren, 'Type','Axes');
% set(ax, 'YTick', [fix(ylim_max(1)/500)*500:500:fix(ylim_max(2)/500)*500], 'YLim', ylim_max)
for i = 1:num_levels
    s.AxesProperties(i).YLimits = ylim_max;
end
xlabel('Time (ms)')
set(gca,'FontSize',10)