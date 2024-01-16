% function plot_averageABR_overlaysingletraces(varargin)
% Plot average ABR waveform (+/- pre and post noise exposure) at a single
% stimulus and frequency in CSV file, with overlaid single traces in gray. 
%
% Default is for user to select file in dialog box when no input parameters
% are specified.
%
% Optional input syntax for specifying input file: plot_averageABR_overlaysingletraces(path, filename)
%
% 9/21/2022 George Liu
% Dependencies: import_averageABRcsv.m, savefig_multiformat.m,
% get_wave1_averageABR.m, get_roc_innerprod_ABR,
% get_thresh_averageABR_liberman.m, get_thresh_averageABR_oghalai.m,
% import_singleABR_frequency.m

opengl('save', 'software') % prevent crashing due to low-level graphics error using Sarah Office computer's graphics card

%% Constants
ABR_PATH = 'D:\George-abr\ABR'; % path to local copy of data, 12-28-21
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures
SAMPLE_PERIOD = 40.96; % sample period of 40.96 us/sample in average trace
SAMPLE_PERIOD_MS = SAMPLE_PERIOD/1000;
COMPRESSION_SINGLEABR_FACTOR = 1;

% Hard coded which stimulus level and frequency to analyze ABRs at
THIS_FREQ = 16000; % Hz
THIS_INTENSITY = 90; % dB

DATE_LABELS = {'Baseline', '24h post', '1 week post'};

% Hard code which mouse to analyze
MOUSE_LABEL = 'b9m8'; % 97 dB noise exposed mouse
MOUSE_DATES = {'20210916', '20210918', '20210925'}; % baseline, post-24 hour, and post-1week time points
ISNOISE = '97 dB'; 

% example of control mouse littermate of 97 dB noise exposed mouse above
MOUSE_LABEL2 = 'b9m9'; % use lowercase
MOUSE_DATES2 = {'20210916', '20210920', '20210926'}; % baseline, post-24 hour, and post-1week time points
ISNOISE2 = 'control'; 

% % example of 99 dB noise exposed mouse
MOUSE_LABEL3 = 'b8m6'; % use lowercase
MOUSE_DATES3 = {'20210911', '20210917', '20210922'}; % baseline, post-24 hour, and post-1week time points
ISNOISE3 = '99 dB';

% Control mouse
MOUSE_LABEL4 = 'b9m10'; % use lowercase
MOUSE_DATES4 = {'20210916', '20210920', '20210925'}; % baseline, post-24 hour, and post-1week time points
ISNOISE4 = 'control';

% Control mouse
MOUSE_LABEL5 = 'b4m4'; % use lowercase
MOUSE_DATES5 = {'20210722', '20210728', '20210803'}; % baseline, post-24 hour, and post-1week time points
ISNOISE5 = 'control';

% 96 dB noise mouse
MOUSE_LABEL6 = 'b6m5'; % use lowercase
MOUSE_DATES6 = {'20210819', '20210822', '20210828'}; % baseline, post-24 hour, and post-1week time points
ISNOISE6 = '96 dB';

% control
MOUSE_LABEL7 = 'b6m2'; % use lowercase
MOUSE_DATES7 = {'20210819', '20210823', '20210828'}; % baseline, post-24 hour, and post-1week time points
ISNOISE7 = 'control';

all_mice_labels = {MOUSE_LABEL, MOUSE_LABEL2, MOUSE_LABEL3, MOUSE_LABEL4, MOUSE_LABEL5, MOUSE_LABEL6, MOUSE_LABEL7};
all_mouse_dates = [MOUSE_DATES; MOUSE_DATES2; MOUSE_DATES3; MOUSE_DATES4; MOUSE_DATES5; MOUSE_DATES6; MOUSE_DATES7];
all_noise_labels = {ISNOISE, ISNOISE2, ISNOISE3, ISNOISE4, ISNOISE5, ISNOISE6, ISNOISE7};

all_mice_labels_2 = {MOUSE_LABEL2, MOUSE_LABEL4};
all_mouse_dates_2 = [MOUSE_DATES2; MOUSE_DATES4];
all_noise_labels_2 = {ISNOISE2, ISNOISE4};

all_mice_labels_3 = {MOUSE_LABEL2};
all_mouse_dates_3 = [MOUSE_DATES2];
all_noise_labels_3 = {ISNOISE2};

N_SINGLETRACES = 40; % number of single traces to plot. Cannot plot all single traces or else MATLAB crashes.

%% Load average trace data
% Iterate over each experiment date

% Show overlayed single traces for 2 mice
these_mice_labels = all_mice_labels_2;
these_mice_dates = all_mouse_dates_2;
these_noise_labels = all_noise_labels_2;

% % Show overlayed single traces for 3 mice
% these_mice_labels = all_mice_labels_3;
% these_mice_dates = all_mouse_dates_3;
% these_noise_labels = all_noise_labels_3;

num_mice = length(these_mice_labels);
num_dates_per_mouse = length(DATE_LABELS);
% M_cache = cell(size(these_mice_dates));
% X_cache = cell(size(these_mice_dates));
x_ms = [];
ylim_max = [0, 0];

% Plot ABRs in one tiled layout
fig = figure;
t = tiledlayout(num_mice, num_dates_per_mouse, 'TileIndexing', 'rowmajor');

for j=1:num_mice
    this_mouse_label = these_mice_labels{j};
    
    for i=1:num_dates_per_mouse
        % Collect ABR data from each date/time point
        this_date = these_mice_dates{j, i};
        
        disp(['  Working on file ', num2str((j-1)*num_dates_per_mouse + i), ' out of ', num2str(numel(these_mice_dates)), ': ', this_mouse_label, ', date ', this_date])

        % Get path to ABR data for experimental date
        filename = [this_date, '_', this_mouse_label, '_abr_left.csv'];
        this_path = fullfile(ABR_PATH, this_date, [this_date, '_', this_mouse_label], 'analyze');

        disp(['Opening ', fullfile(this_path, filename)])
        [M, A_csv, freq_csv] = import_averageABRcsv(filename, this_path);
        num_samples = size(M, 2);
        X = SAMPLE_PERIOD_MS * (1:num_samples); % time in ms
        A_levels = unique(A_csv);
        num_levels = length(A_levels);

        is_freq = freq_csv==THIS_FREQ;
        is_intensity = A_csv==THIS_INTENSITY;

        % ABR trace for specified frequency and intensity
        M2 = M(is_freq & is_intensity, :);

%         M_cache{j, i} = M2;
        
        % 9-21-2022: Obtain single ABR traces at the specified frequency
        % and intensity
        input_pathfile = fullfile(this_path, filename);
        [this_Xcsv, x_ms] = import_singleABR_frequency(input_pathfile, THIS_FREQ, THIS_INTENSITY); % SAMPLES x m/2 matrix
        
        % Compress single trace data to save memory for plotting
        x_ms = x_ms(1:COMPRESSION_SINGLEABR_FACTOR:end);
        this_Xcsv = this_Xcsv(1:COMPRESSION_SINGLEABR_FACTOR:end, :);
        
        % Update y limits
        % Obtain y axis scale
        % Check to make sure bounds of plot are within y limits
        ymax_data = max(this_Xcsv(:));
        ymin_data = min(this_Xcsv(:));
        if ymax_data > ylim_max(2)
            ylim_max(2) = ymax_data;
        end
        if ymin_data < ylim_max(1)
            ylim_max(1) = ymin_data;
        end
               
        % 9-21-2022 Plot single trace ABR waveforms as transparent gray
        % lines
        nexttile(t)
        hold on
        
        m_traces = size(this_Xcsv, 2);
        num_plot_traces = m_traces;
%         num_plot_traces = N_SINGLETRACES;
        cmap = spring(num_plot_traces);
        for k=1:num_plot_traces
            plot(x_ms, this_Xcsv(:, k), 'Color', [cmap(k, 1), cmap(k, 2), cmap(k, 3), 0.8]) % transparency set to 4th color value
        end
        
        % Plot average ABR trace
        y = M2;
        plot(X, y, 'blue', 'LineWidth', 3)
        grid on
        ylim(ylim_max)
        
        % Mark wave 1 peak and following trough
        % Get wave 1 peak and following trough at highest stimulus
        % level
        [peak_pt, trough_pt, amp, lat_peak, lat_trough] = get_wave1_averageABR(X, y);
        CIRCLE_SIZE = 54; % default 36
        scatter(peak_pt(1), peak_pt(2), CIRCLE_SIZE, 'green', 'filled') % peak
        scatter(trough_pt(1), trough_pt(2), CIRCLE_SIZE, 'red', 'filled') % trough
        
        hold off

        % Remove extraneous axis tick marks and x-axis from all but bottom
        % tile
        set(gca,'box','off')
        if j~=num_mice
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'XColor','none')
        end

        % Remove y-axis labels from all but first column
        if i~=1
            set(gca,'yticklabel',[])
        end

        % Set text size of current axis
        set(gca,'FontSize',24)

        % Show frequency as title above top-most tile of each column
%         title([DATE_LABELS{i}, ', ', all_mouse_dates{j, i}])
        if j==1
            title(DATE_LABELS{i})
        end

        % show ylabels for first column only
        if i==1
%             this_ylabel = [all_mice_labels{j}, ' ', all_noise_labels{j}];
            this_ylabel = these_mice_labels{j};
            ylabel(this_ylabel, 'FontSize', 24)
            % rotate y label to make horizontal
            hYLabel = get(gca,'YLabel');
            set(hYLabel,'rotation',90,'VerticalAlignment','bottom', 'HorizontalAlignment','center', 'FontSize', 24)
        end

    end
end
    
% 
% for j=1:num_mice
%     for i=1:num_dates_per_mouse
%         nexttile(t)
%         y = M_cache{j, i};
%         plot(X, y, 'LineWidth', 3)
%         grid on
%         ylim(ylim_max)
% 
%     %     % show ylabels for first column only
%     %     if i==1
%     %         ylabel(newYlabels{i}, 'FontSize', 24)
%     %         % rotate y label to make horizontal
%     %         hYLabel = get(gca,'YLabel');
%     %         set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right', 'FontSize', 24)
%     %     end
% 
%         
%     end
% end


t.TileSpacing = 'none';
% t.TileSpacing = 'tight';
t.Padding = 'tight';
plot_title = ['ABR @ ', num2str(THIS_FREQ), ' Hz, ' num2str(THIS_INTENSITY), ' dB, ' num2str(num_plot_traces), ' single traces'];
title(t, plot_title, 'FontSize', 24)
ylabel(t, 'Amplitude (nV)', 'FontSize', 24)
xlabel(t, 'Time (ms)', 'FontSize', 24)

% Maximize figure window size
fig.WindowState = 'maximized';

% Save figure
%     disp('Saving figure')
% save_file = [this_mouse_label, '_overlaysingletraces'];
num_mice = numel(these_mice_labels);
save_file = '';
for i=1:num_mice
    if i==1
        save_file = these_mice_labels{i};
    else
        save_file = [save_file, '_', these_mice_labels{i}];
    end
end
save_file = [save_file, '_overlaysingletraces_compression', num2str(COMPRESSION_SINGLEABR_FACTOR), 'n', num2str(num_plot_traces)];
savefig_multiformat(gcf, SAVE_PATH, save_file)

disp('Done')

% end