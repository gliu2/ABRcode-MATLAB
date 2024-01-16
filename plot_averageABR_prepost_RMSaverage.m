% function plot_averageABR_prepost_RMSaverage(varargin)
% Correlate wave 1 amplitude variability in control moise with
% characteristics of ABR trace in absence of stimulus. 
%
% Default is for user to select file in dialog box when no input parameters
% are specified.
%
% Optional input syntax for specifying input file: plot_averageABR_RMSsingletraces(path, filename)
%
% Adapted from plot_averageABR_prepost_RMSsingletraces.m
%
% 10/12/2022 George Liu
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
SAMPLING_RATE = 195.3; % samples / msec. Sampling rate of single trial ABR CSV file. Of note, average CSV files have different sampling rate b/c they compress number of samples 1952 -> 244.
SAMPLING_RATE_HZ = SAMPLING_RATE*1000;

% Hard coded which stimulus level and frequency to analyze ABRs at
THIS_FREQ = 16000; % Hz
THIS_INTENSITY = 90; % dB
N_STD_CUTOFF = 5; % number of standard deviations above mean for RMS cutoff for artifact rejection
NUM_ITERATIONS = 6;

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

all_mice_labels_4 = {MOUSE_LABEL2, MOUSE_LABEL4, MOUSE_LABEL5, MOUSE_LABEL7};
all_mouse_dates_4 = [MOUSE_DATES2; MOUSE_DATES4; MOUSE_DATES5; MOUSE_DATES7];
all_noise_labels_4 = {ISNOISE2, ISNOISE4, ISNOISE5, ISNOISE7};

XLIM_RMS = [0, 9000];
YLIM_RMS = [0, 500];

HISTOGRAM_BINS = 50;
N_SINGLETRACES = 40;
WAVEPEAK_SCATTER_CIRCLE_SIZE = 54;
SCATTER_CIRCLE_SIZE = 14; % default 36
FONTSIZE = 14;

SHOW_SCATTER_LEGEND = 'on';
SAVE_FIGURES = 1;


%% Load average trace data
% Iterate over each experiment date

% % Show overlayed single traces for all mice
% these_mice_labels = all_mice_labels;
% these_mice_dates = all_mouse_dates;
% these_noise_labels = all_noise_labels;

% Show overlayed single traces for 2 mice
these_mice_labels = all_mice_labels_2;
these_mice_dates = all_mouse_dates_2;
these_noise_labels = all_noise_labels_2;

% % Show overlayed single traces for 1 mouse
% these_mice_labels = all_mice_labels_3;
% these_mice_dates = all_mouse_dates_3;
% these_noise_labels = all_noise_labels_3;

% % Show overlayed single traces for 1 mouse
% these_mice_labels = all_mice_labels_4;
% these_mice_dates = all_mouse_dates_4;
% these_noise_labels = all_noise_labels_4;

% % Show overlayed single traces for 1 mouse @ 1 date (for evaluation of
% % different cutoffs for artifact rejection and the effect on single traces)
% these_mice_labels = all_mice_labels_3;
% these_mice_dates = all_mouse_dates_3(1);
% these_noise_labels = all_noise_labels_3;

num_mice = length(these_mice_labels);
num_dates_per_mouse = size(these_mice_dates, 2);
% M_cache = cell(size(these_mice_dates));
% X_cache = cell(size(these_mice_dates));
x_ms = [];
ylim_max = [0, 0];

% Colors
mice_colors_unique = hsv(num_mice);

% initialize cache for storing averaged ABR waveforms after artifact
% rejection based on threshold of RMS criterion 
averageABR = cell(num_mice, num_dates_per_mouse);
averageABR_artifactrejection = cell(num_mice, num_dates_per_mouse);
averageABR_artifact_filtered = cell(num_mice, num_dates_per_mouse);
averageABR_noise = cell(num_mice, num_dates_per_mouse);

% Initialize variables for plotting correlation of noise ABR RMS with wave
% 1 amp in stimulus case
wave1amp = zeros(num_mice, num_dates_per_mouse);

% Load average ABR data
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

        averageABR{j, i} = M2; % store average ABR trace
        
        % 10-11-2022: Obtain average ABR in noise at lowest stimlus level
        lowest_stimulus_level = min(A_csv);
        is_lowest_stimulus = A_csv==lowest_stimulus_level;
        averageABR_noise{j, i} = M(is_freq & is_lowest_stimulus, :);
        
        % Obtain y axis scale for plotting average ABR trace later.
        % Check to make sure bounds of plot are within y limits.
        ymax_data = max(M(:));
        ymin_data = min(M(:));
        if ymax_data > ylim_max(2)
            ylim_max(2) = ymax_data;
        end
        if ymin_data < ylim_max(1)
            ylim_max(1) = ymin_data;
        end
        
        

    end
end


%% Plot average ABR waveforms after artifact rejection
% Plot ABRs in one tiled layout
fig = figure;
t = tiledlayout(num_mice, num_dates_per_mouse, 'TileIndexing', 'rowmajor');
for j=1:num_mice
    for i=1:num_dates_per_mouse
        nexttile(t)
        y = averageABR{j, i}; 
        p1 = plot(X, y, 'LineWidth', 3, 'Color', mice_colors_unique(j, :));
        grid on
        ylim(ylim_max)

    %     % show ylabels for first column only
    %     if i==1
    %         ylabel(newYlabels{i}, 'FontSize', 24)
    %         % rotate y label to make horizontal
    %         hYLabel = get(gca,'YLabel');
    %         set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right', 'FontSize', 24)
    %     end

        % Mark wave 1 peak and following trough
        % Get wave 1 peak and following trough at highest stimulus
        % level
        [peak_pt, trough_pt, amp, lat_peak, lat_trough] = get_wave1_averageABR(X, y);

        % 10-11-22: cache wave 1 amplitude to assess for correlation with
        % noise average ABR RMS
        wave1amp(j, i) = amp;
        
        % Show wave 1 peak and trough for original unprocessed average ABR
        % trace
        hold on
        scatter(peak_pt(1), peak_pt(2), WAVEPEAK_SCATTER_CIRCLE_SIZE , 'green', 'filled') % peak
        scatter(trough_pt(1), trough_pt(2), WAVEPEAK_SCATTER_CIRCLE_SIZE , 'red', 'filled') % trough
        hold off
        
        if i==1 && j==1
            legend([p1], {'Raw average'}, 'location', 'southeast', 'FontSize',16)
        end
            
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


t.TileSpacing = 'none';
% t.TileSpacing = 'tight';
t.Padding = 'tight';
plot_title = ['ABR @ ', num2str(THIS_FREQ), ' Hz, ' num2str(THIS_INTENSITY), ' dB'];
title(t, plot_title, 'FontSize', 24)
ylabel(t, 'Amplitude (nV)', 'FontSize', 24)
xlabel(t, 'Time (ms)', 'FontSize', 24)

% Maximize figure window size
fig.WindowState = 'maximized';

% Save figure
%     disp('Saving figure')
base_file = '';
for i=1:num_mice
    if i==1
        base_file = these_mice_labels{i};
    else
        base_file = [base_file, '_', these_mice_labels{i}];
    end
end
if SAVE_FIGURES
    save_file2 = [base_file, '_averageABR_prepost'];
    savefig_multiformat(gcf, SAVE_PATH, save_file2)
end


%% 10-11-2022 Assess for correlation between noise average ABR RMS and wave 1 amplitude
num_data_points = numel(averageABR);
mice_colors_unique = hsv(num_mice);
% Collect average ABR RMS with stimulus
rms_average_abr = zeros(num_mice, num_dates_per_mouse);
rms_no_stimulus_averageABR = zeros(num_mice, num_dates_per_mouse);
max_amp_noise_abr = zeros(num_mice, num_dates_per_mouse);
group_mouse = cell(num_data_points, 1);
group_time = cell(num_data_points, 1);
mice_colors = zeros(num_mice*num_dates_per_mouse, 3);
for j=1:num_mice
    for i=1:num_dates_per_mouse
        ind = i + (j-1)*num_dates_per_mouse;
        rms_average_abr(j, i) = rms(averageABR{j, i});
        group_mouse{ind, 1} = [these_mice_labels{j}, '-', these_noise_labels{j}];
        group_time{ind, 1} = DATE_LABELS{i};
        rms_no_stimulus_averageABR(j, i) = rms(averageABR_noise{j, i});
        max_amp_noise_abr(j, i) = max(averageABR_noise{j, i}) - min(averageABR_noise{j, i});
        
        % add color label
        mice_colors(ind, :) = mice_colors_unique(j, :);
    end
end
group_name = {group_mouse, group_time};

grouped_rms_no_stimulus_averageABR = reshape(rms_no_stimulus_averageABR', [num_data_points, 1]);
grouped_rms_average_abr = reshape(rms_average_abr', [num_data_points, 1]);
grouped_wave1amp = reshape(wave1amp', [num_data_points, 1]);
grouped_max_amp_noise_abr = reshape(max_amp_noise_abr', [num_data_points, 1]);

% Scatter plot colors
mice_symbols = '.*o';

% Scatter plot of wave 1 amp vs noise RMS
fig3 = figure;
t4 = tiledlayout(2,2);
nexttile(t4)
gscatter(grouped_rms_no_stimulus_averageABR, grouped_wave1amp, group_name, mice_colors, mice_symbols, SCATTER_CIRCLE_SIZE, SHOW_SCATTER_LEGEND);
xlabel('RMS of average ABR without stimulus (nV)')
ylabel(['Wave 1 amplitude @ ', num2str(THIS_INTENSITY), ' dB (nV)'])
title(['Wave 1 amplitude vs noise RMS @ ', num2str(THIS_FREQ), ' Hz'])
set(gca,'FontSize', FONTSIZE) % Set text size of current axis

% Scatter plot of average ABR RMS with and without stimulus
nexttile(t4)
gscatter(grouped_rms_no_stimulus_averageABR, grouped_rms_average_abr, group_name, mice_colors, mice_symbols, SCATTER_CIRCLE_SIZE, SHOW_SCATTER_LEGEND);
xlabel('RMS of average ABR without stimulus (nV)')
ylabel(['RMS of average ABR @ ', num2str(THIS_INTENSITY), ' dB (nV)'])
title(['RMS with and without stimulus @ ', num2str(THIS_FREQ), ' Hz'])
set(gca,'FontSize', FONTSIZE) % Set text size of current axis

% Scatter plot of wave 1 amp vs noise max amp
nexttile(t4)
gscatter(grouped_max_amp_noise_abr, grouped_wave1amp, group_name, mice_colors, mice_symbols, SCATTER_CIRCLE_SIZE, SHOW_SCATTER_LEGEND);
xlabel('Max amplitude of average ABR without stimulus (nV)')
ylabel(['Wave 1 amplitude @ ', num2str(THIS_INTENSITY), ' dB (nV)'])
title(['Wave 1 amplitude w/wo stimulus @ ', num2str(THIS_FREQ), ' Hz'])
set(gca,'FontSize', FONTSIZE) % Set text size of current axis

% Scatter plot of average ABR RMS with stimulus vs noise max amp
nexttile(t4)
gscatter(grouped_max_amp_noise_abr, grouped_rms_average_abr, group_name, mice_colors, mice_symbols, SCATTER_CIRCLE_SIZE, SHOW_SCATTER_LEGEND);
xlabel('Max amplitude of average ABR without stimulus (nV)')
ylabel(['RMS of average ABR @ ', num2str(THIS_INTENSITY), ' dB (nV)'])
title(['RMS w stimulus vs peak-peak amp wo @ ', num2str(THIS_FREQ), ' Hz'])
set(gca,'FontSize', FONTSIZE) % Set text size of current axis

t4.TileSpacing = 'tight';
t4.Padding = 'tight';

if SAVE_FIGURES
    save_file3 = [base_file, '_correlate_RMS_wave1amp_withnoise'];
    savefig_multiformat(gcf, SAVE_PATH, save_file3)
end

% 
% fig4 = figure;
% t4 = tiledlayout(num_mice, num_dates_per_mouse, 'TileIndexing', 'rowmajor');
% 
% % Determine max RMS for y-limits
% max_rms = cellfun(@max, singletraces_rms_cache);
% % max_rms_all = max(max_rms(:));
% max_rms_all = 8000;
% for j=1:num_mice
%     for i=1:num_dates_per_mouse
%         % Plot RMS vs single trace index
%         nexttile(t)
%         y_rms = singletraces_rms_cache{j, i}; 

disp('Done')

% end