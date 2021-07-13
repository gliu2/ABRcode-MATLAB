% Script to plot and analyze DPOAE data 
%
% 7/12/2021 George Liu
% Dependencies: import_DPOAEcsv.m


%% Load data
[filename,path] = uigetfile('*.csv');
[M, A_csv, freq_csv, f1, f2, sample_period] = import_DPOAEcsv(filename, path);
num_samples = size(M, 2);
num_levels = size(M, 1);
plot_title = ['DPOAE @ ', num2str(freq_csv(1)), ' Hz'];

% X = sample_period / 1000 * (1:num_samples); 
X = 97656.3/2048* (1:num_samples);
dB = cell(num_levels, 1);
dB(:) = {'dB'};
[A_descending, I] = sort(A_csv, 'descend');
M_sorted = M(I, :);
A_descending_cell = cellstr(num2str(A_descending));
newYlabels=cellfun(@(x,y) [x  ' ' y], A_descending_cell, dB, 'un', 0);

%% Plot DPOAE measurements in time domain in vertical stack
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
    s.AxesProperties(i).XLimits = [0,10000];
end
xlabel('Frequency (Hz)')
set(gca,'FontSize',10)