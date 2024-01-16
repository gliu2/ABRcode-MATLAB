% function plot_CAP_amp_change
% Plot data in Excel sheet of metrics (CAP amplitude, CM amplitude,
% normalized CAP/CM amplitude) vs timepoint for each mouse (row in Excel sheet).
%
% 8/22/2023 George Liu
% Dependencies: savefig_multiformat.m


SAVE_PATH = 'd:\users\admin\Documents\George\Results_analyzed'; % path for saving figures

disp('Select file to open')
[filename,path] = uigetfile('*.xlsx', 'Select data file');
disp(['Opening ', fullfile(path, filename)])

T = readtable(fullfile(path, filename));

%%
data_cols = {[5, 6, 7], [8, 9, 10], [11, 12, 13]};
data_labels = ["CM amplitude", "CAP wave 1 amplitude", "CAP wave 1 amp / CM amp"];
% time_labels = ["Baseline", "Post 10 min", "Post 60 min"];
timepoints = [-40, 10, 60];
n_sets = length(data_cols);
yunits = ["CM amplitude (nV)", "CAP wave 1 amplitude (nV)", "Normalized amplitude (a.u.)"];
num_rows = size(T, 1);
c = parula(num_rows);

str1 = datestr(T{:, 1});
str2 = T{:, 2};
num_rows = size(T, 1);
leg_str = cell(num_rows, 1);
for i=1:num_rows
    leg_str{i} = [str1(i, :), str2{i}];
end

for i = 1:n_sets
    fig = figure;
    hold on
    for j = 1:size(T, 1)
        plot(timepoints, T{j, data_cols{i}}, 'LineWidth', 3, 'Color', c(j, :))
    end
    
    % Set text size of current axis
    set(gca,'FontSize',24)
    
    ylabel(yunits(i), 'FontSize', 24)
    xlabel("Time (min)", 'FontSize', 24)
    title(data_labels(i), 'FontSize', 24)
    xlim([min(timepoints), max(timepoints)])
    
    legend(leg_str, 'Location', 'eastoutside')
    
    % Maximize figure window size
    fig.WindowState = 'maximized';  
    
    % Save figure
    save_file = convertStringsToChars(data_labels(i));
    save_file = regexprep(save_file, ' ', '_');
    save_file = regexprep(save_file, '_/_', '_');
    savefig_multiformat(gcf, SAVE_PATH, [save_file, '_CAPCMchange'])
end

disp('Done')