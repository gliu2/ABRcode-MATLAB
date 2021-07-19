% Script to plot average ABRs in CSV file in stacked plot
%
% 7/11/2021 George Liu
% Dependencies: import_averageABRcsv.m, savefig_multiformat.m

%% Constants
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures
SAMPLE_PERIOD = 40.96; % sample period of 40.96 us/sample in average trace

%% Load average trace data
[filename,path] = uigetfile('*.csv');
disp(['Opening ', fullfile(path, filename)])
[M, A_csv, freq_csv] = import_averageABRcsv(filename, path);
num_samples = size(M, 2);
X = SAMPLE_PERIOD / 1000 * (1:num_samples); 

% Analyze DPOAEs for each frequency separately
[freq_unique, ~, ic] = unique(freq_csv);
n_freq = length(freq_unique);

for ff = 1:n_freq
    this_freq = freq_unique(ff);
    disp(['Working on ', num2str(this_freq), ' Hz...'])
    
    % amplitudes and DPOAE values for this frequency only
    A_csv2 = A_csv(ic==ff, :);
    M2 = M(ic==ff, :);
    
    num_levels = size(M2, 1);
    dB = cell(num_levels, 1);
    dB(:) = {'dB (nV)'};

    [A_descending, I] = sort(A_csv2, 'descend');
    M_sorted = M2(I, :);
    A_descending_cell = cellstr(num2str(A_descending));
    newYlabels=cellfun(@(x,y) [x  ' ' y], A_descending_cell, dB, 'un', 0);

    %% Plot ABRs in vertical stack
    figure
    plot_title = ['ABR @ ', num2str(this_freq), ' Hz'];
    s = stackedplot(X, M_sorted', 'Title', plot_title, 'DisplayLabels', newYlabels, 'LineWidth', 1.5);
    % match y axis scale for all plots
%     ylim_max = [0, 0];
%     for j = 1:num_levels
%         this_ylim_max = s.AxesProperties(j).YLimits;
%         ylim_max(1) = min([ylim_max(1), this_ylim_max(1)]);
%         ylim_max(2) = max([ylim_max(2), this_ylim_max(2)]);
%     end
    ylim_max = [-1200, 1200];
    % CHeck to make sure bounds of plot are within y limits
    ymax_data = max(M(:));
    ymin_data = min(M(:));
    if ymax_data > ylim_max(2)
        ylim_max(2) = ymax_data;
    end
    if ymin_data < ylim_max(1)
        ylim_max(1) = ymin_data;
    end
    
    % ax = findobj(s.NodeChildren, 'Type','Axes');
    % set(ax, 'YTick', [fix(ylim_max(1)/500)*500:500:fix(ylim_max(2)/500)*500], 'YLim', ylim_max)
    for i = 1:num_levels
        s.AxesProperties(i).YLimits = ylim_max;
    end
    xlabel('Time (ms)')
    set(gca,'FontSize',10)
    
    % Save figure
    disp('Saving figure')
    [~, save_file, ~] = fileparts(filename);
    save_file_freq = [save_file, num2str(this_freq), 'hz'];
    savefig_multiformat(gcf, SAVE_PATH, save_file_freq)
    
end

disp('Done')