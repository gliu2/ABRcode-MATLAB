% function thresholdpaper_fig2(varargin)
% Plot estimated threshold vs manually calculated threshold
%
% 8/27/2023 George Liu
% Dependencies: same_xaxes.m, same_yaxes.m

%% Constants
LOAD_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures
SAVE_PATH = 'd:\users\admin\Documents\George\Results_analyzed'; % path for saving figures

SAVE_FIGURES = 1;

NUM_CHANNELS = 2; % 1 = ABR, 2 = CAP
CHANNEL_KEY = {'ABR', 'CAP'};
NUM_UNDERSCORES = 4; % Expected number of underscores in filename of .mat output results file
NUM_DASHES = 0;
n_freq = 3;

% Pre-specified list of control data to analyze
CONTROL_TUBE_FILES = { ...
    '20230725_control1'; ...
    '20230726_control1'; ...
    '20230726_control2' ...
    };

CONTROL_AP_FILES = { ...
    '20230727_control1'; ...
    '20230727_control2'; ...
    '20230727_control3'; ...
    '20230813_control2' ...
    };

CONTROL_CNQX_FILES = { ...
    '20230801_cnqx1'; ...
    '20230801_cnqx2'; ...
    '20230802_cnqx1' ...
    };

EXPERIMENT_GROUP2 = [CONTROL_TUBE_FILES; CONTROL_AP_FILES; CONTROL_CNQX_FILES];

%% Obtain file metadata of result '.mat' files in Results folder.
listing = dir(fullfile(LOAD_PATH, '*.mat'));
allnames_all = {listing.name}'; % column cell array of filenames

underscore_loc_all = strfind(allnames_all, '_');
is_enough_underscores = cellfun(@numel, underscore_loc_all) >= NUM_UNDERSCORES;
dash_loc_all = strfind(allnames_all, '-');
is_no_dashes = cellfun(@numel, dash_loc_all) == NUM_DASHES;
is_results = is_enough_underscores & is_no_dashes;

allnames = allnames_all(is_results);
underscore_loc = underscore_loc_all(is_results);
n_files = length(allnames);
dates = cell(n_files, 1);
label = cell(n_files, 1);
timepoint = cell(n_files, 1);
channel = cell(n_files, 1);

for i = 1:n_files
    dates(i) = extractBetween(allnames{i}, 1, underscore_loc{i}(1) - 1); 
    label(i) = extractBetween(allnames{i}, underscore_loc{i}(1) + 1, underscore_loc{i}(2) - 1); 
    timepoint(i) = extractBetween(allnames{i}, underscore_loc{i}(2) + 1, underscore_loc{i}(3) - 1); 
    channel(i) = extractBetween(allnames{i}, underscore_loc{i}(3) + 1, underscore_loc{i}(4) - 1); 
end

%% Extract results 
% select dataset
this_filenames = EXPERIMENT_GROUP2;

% Get selected filename dates and labels
n_files_analyze = length(this_filenames);
analyze_underscore_loc = strfind(this_filenames, '_');

this_date = cell(n_files_analyze, 1);
this_label = cell(n_files_analyze, 1);
is_date_and_label = zeros(n_files, 1);
for i = 1:n_files_analyze
    this_date{i} = extractBefore(this_filenames{i}, analyze_underscore_loc{i}(1)); 
    this_label{i} = extractAfter(this_filenames{i}, analyze_underscore_loc{i}(1));  

    is_date = strcmp(dates, this_date{i});
    is_label = strcmp(label, this_label{i});

    is_date_and_label = is_date_and_label | (is_date & is_label);
end

% Construct table with extracted results
% theseResults = 
analyze_names = allnames(is_date_and_label);
Date = dates(is_date_and_label);
Mouse = label(is_date_and_label);
Time = timepoint(is_date_and_label);
Channel = channel(is_date_and_label);
n_analyze = length(analyze_names);
Threshold_8khz = zeros(n_analyze , 1);
Threshold_16khz = zeros(n_analyze , 1);
Threshold_32khz = zeros(n_analyze , 1);
ThresholdOgh_8khz = zeros(n_analyze , 1);
ThresholdOgh_16khz = zeros(n_analyze , 1);
ThresholdOgh_32khz = zeros(n_analyze , 1);


for i = 1:n_analyze 
    load(fullfile(LOAD_PATH, analyze_names{i}))

    Threshold_8khz(i) = threshold_D(1);
    Threshold_16khz(i) = threshold_D(2);
    Threshold_32khz(i) = threshold_D(3);

%     if strcmp(Channel{i}, CHANNEL_KEY{2}) % CAP
%         Threshold_8khz(i) = threshold_D(1);
%         Threshold_16khz(i) = threshold_D(2);
%         Threshold_32khz(i) = threshold_D(3);
%     end

    ThresholdOgh_8khz(i) = threshold(1);
    ThresholdOgh_16khz(i) = threshold(2);
    ThresholdOgh_32khz(i) = threshold(3);

%     'y_avg_cache_descending', ...
%         'y_avg_cm_cache_descending', 'all_freqs', 'A_descending', ...
%         'threshold', 'metric', 'threshold_D', 'metric_D', ... % don't save 'X_csv_cache' because too large
%         'amp_cm_8khz90db', ...
%         'wave1peak_cache', ...
%         'wave1trough_cache', ...
%         'wave1amp_cache', ...
%         'wave1lat_peak_cache', ...
%         'wave1lat_trough_cache', ...
%         'dist_wave1amp', ...
%         'wave1std' ...
%         )
end

extractedResults = table(Date, Mouse, Time, Channel, ...
    Threshold_8khz, Threshold_16khz, Threshold_32khz, ...
    ThresholdOgh_8khz, ThresholdOgh_16khz, ThresholdOgh_32khz ...
    );

%% Plot analysis
% for j = 1:NUM_CHANNELS % 1 = ABR, 2 = CAP
for j=1
    % Group average of threshold change
    is_channel = strcmp(extractedResults.Channel, CHANNEL_KEY{j}); 
    is_pre = strcmp(extractedResults.Time, 'pre');

    extractedResults_pre = extractedResults(is_channel & is_pre, :);
    threshold_cap_pre = extractedResults{is_channel & is_pre, ["Threshold_8khz", "Threshold_16khz", "Threshold_32khz"]};
end
   
%% Load data from Excel file (includes human threshold estimates)
PATH_EXCEL_THRESHOLD = 'd:\users\admin\Documents\George\Results_manual\analyze_CAPABR_thresholds_8-30-23.xlsx'; % George measurements
PATH_EXCEL_THRESHOLD_SRI = 'D:\users\admin\Documents\George\CAPABR\analyze_CAPABR_thresholds_template_8-31-23_sri.xlsx';
T = readtable(PATH_EXCEL_THRESHOLD);
T_human = readtable(PATH_EXCEL_THRESHOLD_SRI);
x_limits = [0, 100];
y_limits = [-50, 50];

for i = 1:NUM_CHANNELS
% for i = 2
    rows_for_channel = {11:20, 1:10};
    % Plot Bland-Altman of baseline mouse CAP thresholds - human vs algorithm
    thresh_D = T{rows_for_channel{i}, 5:7}; % 10 x 3 matrix, each column is different frequency
    thresh_human1 = T{rows_for_channel{i}, 11:13}; % 10 x 3 matrix, each column is different frequency
    thresh_human2 = T_human{rows_for_channel{i}, 5:7}; % 10 x 3 matrix, each column is different frequency
%     thresh_human = (thresh_human1 + thresh_human2)/2;
    thresh_D = thresh_human1;
    thresh_human = thresh_human2;

    N_SUBJECTS = 10;
    CIRCLE_SIZE = 36;
    FONT_SIZE = 14;

    fig = figure;
    t = tiledlayout(1, n_freq, 'TileIndexing', 'columnmajor');
    axesHandle = zeros(n_freq, 1);
    mean_dif_cache = zeros(n_freq, 1);
    std_dif_cache = zeros(n_freq, 1);
    for j=1:n_freq
        axesHandle(j) = nexttile(t);
        hold on

        these_thresh = [thresh_D(:, j), thresh_human(:, j)];
        thresh_avg = nanmean(these_thresh, 2);
        thresh_dif = these_thresh(:, 1) - these_thresh(:, 2);

        g1 = 1:N_SUBJECTS;
        c = hsv(N_SUBJECTS);

        gscatter(thresh_avg, thresh_dif, g1, c, '.', CIRCLE_SIZE, 'doleg', 'off')

        xlabel('Mean threshold (dB SPL)')
%         ylabel('Difference (algorithm - human) (dB)')
        ylabel('Difference (dB)')
        title(sprintf('%.0f kHz', all_freqs(j)/1000))
        set(gca, 'FontSize', FONT_SIZE)

        mean_dif = nanmean(thresh_dif);
        std_dif = nanstd(thresh_dif);
        yline(mean_dif, '-', sprintf('MEAN: %.0f dB', mean_dif), 'HandleVisibility','off', 'LineWidth', 2, 'FontSize', FONT_SIZE)
        yline(mean_dif + std_dif, '--', sprintf('+SD: %.0f dB', mean_dif + std_dif), 'HandleVisibility','off', 'LineWidth', 2,  'FontSize', FONT_SIZE)
        yline(mean_dif - std_dif, '--', sprintf('-SD: %.0f dB', mean_dif - std_dif), 'HandleVisibility','off', 'LineWidth', 2, 'FontSize', FONT_SIZE)

        hold off
        
        xticks([0:10:100])
        yticks([-50:10:50])
        xlim(x_limits)
        ylim(y_limits)
        axis square
%         axis equal

        % Cache mean_dif and std_dif
        mean_dif_cache(j) = mean_dif;
        std_dif_cache(j) = std_dif;
    end
%     title(t, ['Bland-Altman plot of baseline mouse ', CHANNEL_KEY{i}], 'FontSize', 24)
    
    ylabel(t, ['Mouse ', CHANNEL_KEY{i}], 'FontSize', 24)
    t.TileSpacing = 'loose';
    t.Padding = 'tight';

    % Maximize figure window size
    fig.WindowState = 'maximized'; 
    
    same_xaxes(axesHandle)
    same_yaxes(axesHandle)
    
    %% Save figure
    if SAVE_FIGURES
        save_file = 'blandaltman';
%         savefig_multiformat(gcf, SAVE_PATH, [save_file, '_', CHANNEL_KEY{i}, '_thresholdpaper_Fig2'])
%         savefig_multiformat(gcf, SAVE_PATH, [save_file, '_', CHANNEL_KEY{i}, '_thresholdpaper_Fig2_sri'])
%         savefig_multiformat(gcf, SAVE_PATH, [save_file, '_', CHANNEL_KEY{i}, '_thresholdpaper_Fig2_sri_gl'])
        savefig_multiformat(gcf, SAVE_PATH, [save_file, '_', CHANNEL_KEY{i}, '_thresholdpaper_Fig2_humans'])
    end
    
end

%% Plot bar plots
PATH_EXCEL_BARDATA = 'D:\users\admin\Documents\George\CAPABR\analyze_CAPABR_thresholds_blandaltman_barplots_9-1-23.xlsx';
T_bar = readtable(PATH_EXCEL_BARDATA);

COL_MEANDIF_ABR = 3;
COL_MEANDIF_CAP = 4;
COL_STDDIF_ABR = 5;
COL_STDDIF_CAP = 6;

COL_MEANDIF_CHANNEL = [COL_MEANDIF_ABR, COL_MEANDIF_CAP];
COL_STDDIF_CHANNEL = [COL_STDDIF_ABR, COL_STDDIF_CAP];

N_MICE = 10;
for i = 1:NUM_CHANNELS
% for i = 2
    mean_dif_D = T_bar{1:3, COL_MEANDIF_CHANNEL(i)}; 
    mean_dif_human = T_bar{4:6, COL_MEANDIF_CHANNEL(i)};
    stddif_D = T_bar{1:3, COL_STDDIF_CHANNEL(i)};
    stddif_human = T_bar{4:6, COL_STDDIF_CHANNEL(i)};
    
    mean_dif_D_all = mean(mean_dif_D);
    mean_dif_human_all = mean(mean_dif_human);
    % Pool standard deviations across all 3 frequencies
%     stddif_D_all = mean(stddif_D);
    [n_2samp, mean_dif_D_2samp, stddif_D_2samp] = pooledmeanstd(N_MICE, mean_dif_D(1), stddif_D(1), N_MICE, mean_dif_D(2), stddif_D(2));
    [~, ~, stddif_D_all] = pooledmeanstd(n_2samp, mean_dif_D_2samp, stddif_D_2samp, N_MICE, mean_dif_D(3), stddif_D(3));
%     stddif_human_all = mean(stddif_human);
    [n_2samp, mean_dif_human_2samp, stddif_human_2samp] = pooledmeanstd(N_MICE, mean_dif_human(1), stddif_human(1), N_MICE, mean_dif_human(2), stddif_human(2));
    [~, ~, stddif_human_all] = pooledmeanstd(n_2samp, mean_dif_human_2samp, stddif_human_2samp, N_MICE, mean_dif_human(3), stddif_human(3));
    
    mean_dif_D_extended = [mean_dif_D; mean_dif_D_all];
    mean_dif_human_extended = [mean_dif_human; mean_dif_human_all];
    stddif_D_extended = [stddif_D; stddif_D_all];
    stddif_human_extended = [stddif_human; stddif_human_all];
    
    mean_dif_combined = [mean_dif_D_extended, mean_dif_human_extended];
    stddif_combined = [stddif_D_extended, stddif_human_extended];
    
    
    figure
    b = bar(abs(mean_dif_combined), 'grouped');
    hold on
    % Calculate the number of groups and number of bars in each group
    [ngroups,nbars] = size(mean_dif_combined);
    % Get the x coordinate of the bars
    x = nan(nbars, ngroups);
    for j = 1:nbars
        x(j,:) = b(j).XEndPoints;
    end
    % Plot the errorbars
    errorbar(x', abs(mean_dif_combined), zeros(size(stddif_combined)), stddif_combined, 'k', 'linestyle','none');
    hold off
    
    set(gca, 'FontSize', FONT_SIZE)
    ylabel('Mean threshold difference (dB)')
    xlabel('Frequency (kHz)')
    title(CHANNEL_KEY{i})
    hLg=legend({'Algorithm vs humans', 'Human vs human'},'Location','northeast');
    group_names = {'8', '16', '32', 'All'};
    set(gca,'xticklabel', group_names)
    
    ylim([0, 25])
     
    % Save figure
    if SAVE_FIGURES
        save_file = 'blandaltman';
        savefig_multiformat(gcf, SAVE_PATH, [save_file, '_', CHANNEL_KEY{i}, '_thresholdpaper_Fig2c_bars'])
    end
    
   
end
