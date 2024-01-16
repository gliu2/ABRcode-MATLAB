% function analyze_CAPABR_results(varargin)
% Analyze CAP and ABR .mat results output by plotstack_average_CAPABR.m.
%
% 8/27/2023 George Liu
% Dependencies: ***

%% Constants
LOAD_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures
SAVE_PATH = 'd:\users\admin\Documents\George\Results_analyzed'; % path for saving figures

SAVE_FIGURES = 1;

NUM_CHANNELS = 2; % 1 = ABR, 2 = CAP
CHANNEL_KEY = {'ABR', 'CAP'};
NUM_UNDERSCORES = 4; % Expected number of underscores in filename of .mat output results file
NUM_DASHES = 0;

% Pre-specified list of control data to analyze
CONTROL_NOTUBE_FILES = { ...
    '20230823_control1'; ...
%     '20230823_control2'; ... % no baseline hearing at 32 kHz
    '20230823_control3'; ...
    '20230823_control4'; ...
    '20230824_control1' ...
    };

CONTROL_TUBE_FILES = { ...
    '20230725_control1'; ...
    '20230726_control1'; ... % no baseline hearing at 32 kHz
    '20230726_control2'; ...
    '20230828_tube2' ...
    };

CONTROL_AP_FILES = { ...
    '20230727_control1'; ...
    '20230727_control2'; ...
%     '20230727_control3'; ... % no baseline hearing at 32 kHz
    '20230813_control2'; ...
%     '20230903_control1'; ... % Exclude because waves change shape significantly
    '20230903_control2'; ...
%     '20230903_control3'; ... % Exclude because waves change shape significantly
    };

CNQX_FILES = { ...
%     '20230801_cnqx1'; ... % no baseline hearing at 32 kHz
    '20230801_cnqx2'; ...
    '20230802_cnqx1'; ...
    '20230828_cnqx1'; ...
    '20230828_cnqx3'; ...
    '20230829_cnqx1'; ...
    '20230829_cnqx2'; ...
    '20230831_cnqx1'; ...
    '20230831_cnqx2'; ...
    '20230831_cnqx3' ...
    };

BAD_AP_FILES = { ...
    '20230721_control1'; ...
    '20230722_control1'; ...
    '20230723_control1'; ... 
    '20230724_control1'; ... 
    '20230724_control2' ... 
    };

EXPERIMENT_GROUP = {CONTROL_NOTUBE_FILES, CONTROL_TUBE_FILES, CONTROL_AP_FILES, CNQX_FILES, BAD_AP_FILES};
EXPERIMENT_TITLES = {'No tube', 'Tube', 'Artificial perilymph', '100 uM CNQX', 'Bad AP'};
EXPERIMENT_SAVEFILENAME = {'control_notube', 'control_tube', 'control_ap', 'cnqx', 'bad_ap'};

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
for k = 1:length(EXPERIMENT_GROUP)
    this_filenames = EXPERIMENT_GROUP{k};
    % this_filenames = CONTROL_TUBE_FILES;
    % this_filenames = CONTROL_AP_FILES;
%     this_filenames = CONTROL_CNQX_FILES;

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
    for j = 1:NUM_CHANNELS % 1 = ABR, 2 = CAP
        % Group average of threshold change
        is_channel = strcmp(extractedResults.Channel, CHANNEL_KEY{j}); 
        if strcmp(EXPERIMENT_SAVEFILENAME{k}, 'control_notube')
            is_pre = strcmp(extractedResults.Time, 'pre');
            is_post10 = strcmp(extractedResults.Time, 'pre2');
            is_post30 = strcmp(extractedResults.Time, 'pre3');
            is_post60 = strcmp(extractedResults.Time, 'pre4');            
        else
            is_pre = strcmp(extractedResults.Time, 'pre');
            is_post10 = strcmp(extractedResults.Time, 'post10');
            is_post30 = strcmp(extractedResults.Time, 'post30');
            is_post60 = strcmp(extractedResults.Time, 'post60');
        end

        threshold_cap_pre = extractedResults{is_channel & is_pre, ["Threshold_8khz", "Threshold_16khz", "Threshold_32khz"]};
        threshold_cap_post10 = extractedResults{is_channel & is_post10, ["Threshold_8khz", "Threshold_16khz", "Threshold_32khz"]};
        threshold_cap_post30 = extractedResults{is_channel & is_post30, ["Threshold_8khz", "Threshold_16khz", "Threshold_32khz"]};
        threshold_cap_post60 = extractedResults{is_channel & is_post60, ["Threshold_8khz", "Threshold_16khz", "Threshold_32khz"]};

        % Plot average threshold change
        threshold_cap_timepoints = {threshold_cap_pre, threshold_cap_post10, threshold_cap_post30, threshold_cap_post60};
        timepoint_labels = ["Pre", "10 min", "30 min", "60 min"];
        n_timepoints = length(timepoint_labels);
        colors = [[0.5, 0.5, 0.5]; flipud(copper(n_timepoints-1))];
        x = all_freqs/1000; % kHz
        n_mice = [];

        XLIM = [7, 35];
        YLIM = [0, 100];

        fig = figure;
        hold on
        for i = 1:n_timepoints
            y = threshold_cap_timepoints{i};
            y_mean = nanmean(y, 1);
            y_ste = nanstd(y, 1)./sqrt(sum(~isnan(y),1));

            errorbar(x, y_mean, y_ste, 'LineWidth', 3, 'Color', colors(i, :))

            n_mice = max([n_mice, size(y,1)]); % Obtain number of mice used in experiments. Some mice have invalid 32 kHz hearing results.
        end
        % Plot asterisk above frequency if threshold is significantly
        % different 1 hour after baseline
        P_CRIT = 0.05;
        n_freq = length(all_freqs);
        for i=1:n_freq
            [~, p] = ttest2(threshold_cap_pre(:, i), threshold_cap_post60(:, i), 'Tail', 'both'); % compare control and noise synapse means
            if p < P_CRIT
%                 disp(['Significant difference: ', num2str(x(i)), ' kHz'])
                text(x(i), y_mean(i) + y_ste(i) + 10, '*', 'Color', 'k', 'FontSize', 28, 'FontWeight', 'bold')
            end
        end
        hold off
        legend(timepoint_labels, 'Location', 'southeast') 
        legend boxoff   

        set(gca, 'XScale', 'log');
        xticks(x)
        xlim(XLIM)
        ylim(YLIM)

        xlabel('Frequency (kHz)', 'FontSize', 24)
        ylabel('Threshold (dB SPL)', 'FontSize', 24)

        title([EXPERIMENT_TITLES{k}, ', ', CHANNEL_KEY{j}, ', N=', num2str(n_mice)], 'FontSize', 24, 'Interpreter', 'none')
        % title(['Tube, N=', num2str(size(y,1))], 'FontSize', 24, 'Interpreter', 'none')
        % title(['Artificial perilymph, N=', num2str(size(y,1))], 'FontSize', 24, 'Interpreter', 'none')
%         title(['100 uM CNQX', ', N=', num2str(n_mice)], 'FontSize', 24, 'Interpreter', 'none')

        % Set text size of current axis
        set(gca,'FontSize',24)

        % Maximize figure window size
        fig.WindowState = 'maximized';
        axis square

        %% Save figure after manually resizing window
        if SAVE_FIGURES
            % Save figure
            %     disp('Saving figure')
        %     [~, save_file, ~] = fileparts(filename_prestem);
            save_file = [EXPERIMENT_SAVEFILENAME{k}, '_group'];
        %     save_file = 'control_tube_group';    
        %     save_file = 'control_ap_group_CAP';
%             save_file = 'cnqx_group';
            save_channel = CHANNEL_KEY{j};
            savefig_multiformat(gcf, SAVE_PATH, [save_file, '_', save_channel, '_thresholdshift_9-1-23'])
        end
    end
end