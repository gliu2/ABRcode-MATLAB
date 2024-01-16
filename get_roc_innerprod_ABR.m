function signal_metrics_thresholds = get_roc_innerprod_ABR(varargin)
% function [auc_cache, ip_mean_cache, ip_ste_cache, thresh_cache, thresh_criterion, d_prime_cache_rms, thresh_d_prime_rms, d_prime_cache_z, thresh_d_prime_z] = get_roc_innerprod_ABR(varargin)
% 
% Calculate inner product (mean) and ROC-AUC level functions for single ABR
% traces.
%
% For a selected average ABR CSV file, get path of parent folder
% containing single trace CSV files of same filename prefix. Analyze with
% raw inner product function. 
% 
% Inputs (optional):
%       - varargin{1}: full file path to average ABR CSV, related to single
%       ABR CSV files. 
%       - varargin{2}: logical value (true/false) for whether to plot
%       figures
%       - varargin{3}: optional specification of wave 1 peak and trough latencies
%       to perform single trace ABR analysis focused around wave 1. Size
%       (n_freq x num_levels, 4) matrix.
%           -Column 1: wave 1 peak latency (ms)
%           -Column 2: wave 1 trough latency (ms)
%           -Column 3: stimulus frequency (Hz) 
%           -Column 4: stimulus level in descending order (dB) 
%
% Outputs: 
%   signal_metrics_thresholds - cell encapsulating signal detection sensitivity metrics and thresholds listed
%   below. Cell of size (n_windows, 1), where n_windows is the number of
%   different windows used to analyze subset of trace timepoints (if optional varargin{3} 
%   input is specified). Default is to analyze only the whole trace.
%       - auc_cache:    ROC-AUC values for each set of single ABR traces.
%                       Matrix of size (n_freq, num_levels).
%       - ip_mean_cache: Mean inner product value of all single traces'
%                       inner products with coherent average at same stimulus 
%                       level. Matrix of size (n_freq, num_levels). 
%       - ip_ste_cache: Standard error of inner product mean. Matrix of size (n_freq, num_levels). 
%       - thresh_cache: Estimated ABR thresholds, based on Wilcoxin
%                       sum-rank test p-val < 0.5 for inner product
%                       distribution with noise distribution. Vector of length (n_freq). In order of freqs_unique. 
%       - thresh_criterion: ABR thresholds calculated to determine lowest stimulus
%                       level where mean inner product is greater than
%                       criterion set by mean innerproduct + 3 standard
%                       deviations with no stimulus.
% 
%       12/6/21 - removed these outputs (were present in prior versions of
%       code)
%       - auc_lower_cache:    ROC-AUC lower bound confidence estimate based on bootstrapping.
%                       Matrix of size (n_freq, num_levels).
%       - auc_upper_cache:    ROC-AUC upper bound confidence estimate based on bootstrapping.
%                       Matrix of size (n_freq, num_levels).
%
% Example usage: 
% [auc_cache, ip_mean_cache, thresh_cache, A_desc, freqs_unique] = get_roc_innerprod_ABR;
%
% George Liu
% Last edited 9/21/2022
% 12/8/21 - analyze single trial ABRs using full trace and windowed trace
% in parallel (instead of sequentially) to reduce runtime from loading
% single traces twice. -> Output is cell encapsulating both results.
% 
% 
% 12/6/21 - add d' metric (difference between the means of the signal and noise
% inner product distribution histograms, normalized by the RMS of their
% standard deviations) as another output. This is another candidate metric
% of signal detection sensitivity, in addition to AUC (captures information
% in the mean and width of signal IP distribution histograms unlike AUC and IP mean).
% - removed bootstrapping for AUC upper and lower estimates to speed up
% code
% 
% 11/12/21 - added threshold_criterion output determined using
% criterion level similar to Oghalai method.
%
% 10/22/21 - made time window option input an array of nx2
% rather than 2 element vector, so that time window can be different for
% each stimulus level
% 
% 9/8/21 
% - added 3rd input option to allow innerprod analysis on selected time window of ABR trace.
% - added criteria that ABR threshold require p-val < 0.5 at AND for all
% stimulus levels above threshold.
% - added bootstrapping to estimate AUC lower and upper confidence intervals for standard errors
% - calculate standard error for inner product mean
% 
% 9/1/21 - save NaN values for missing single trace ABR files
% (eg if first set of single trace files were accidentally during saving of average ABR files).
%
% Dependencies: import_ABRcsv.m, get_innerprod_raw_ABR.m, get_innerprod_matrix_ABR.m,
% merge_singletraceABR_polarities.m, xgiveny_avetemplate.m

% [filename, path] = uigetfile('*.csv');
% [X_csv, A_csv, freq_csv] = import_ABRcsv(filename, path);
% this_ip_lone = get_innerprod_raw_ABR(X_csv);
% figure, histogram(this_ip_lone)

%% Constants
MAX_THRESHOLD = 95; % dB; arbitrary max threshold if ABRs are noisy at all stimuli levels
SAMPLING_RATE = 195.3; % samples / msec. Sampling rate of single trial ABR CSV file. Of note, average CSV files have different sampling rate b/c they compress number of samples 1952 -> 244.
TRIM_ABR = false;
N_SINGLE_TRACES = 1026;
TIME_AROUNDPEAK = 0.2; % ms, time around peak and trough, for calculating inner product around wave 1
TIME_BEFOREPEAK = 0.2; % start time window this many ms before wave 1 peak - added 12/6/21
TIME_WINDOW_LENGTH = 1.2; % All time windows are this long (ms) from constant time before wave 1 peak
X_LIM_WHOLETRACE = [-5, 5]*10^9; % maximum inner product value (nV^2) for plotting histograms of inner product distributions using the whole trace (no window). if inner products are not normalized 

%%
% file_pre = '20210823_b6m2_abr_left';
% file_ext = '.csv';
% file_pattern = [file_pre, '*', file_ext];
% selpath = uigetdir();

% Select average ABR CSV file, for which single traces you would like to
% analyze
if nargin==0
    [selfile, selpath_analyze] = uigetfile('*.csv');
    input_pathfile = fullfile(selpath_analyze, selfile);
    PLOT_FIGURES_DEFAULT = true;
elseif nargin==1
    input_pathfile = varargin{1};
    PLOT_FIGURES_DEFAULT = false;
elseif nargin==2
    input_pathfile = varargin{1};
    PLOT_FIGURES_DEFAULT = varargin{2};
elseif nargin==3
    input_pathfile = varargin{1};
    PLOT_FIGURES_DEFAULT = varargin{2};
    wave1_latency = varargin{3}; % input latencies of wave 1 peak and trough (ms)
    TRIM_ABR = true;
else
    disp('Warning! Number of input arguments cannot be more than 3.')
end
    
[selpath_analyze, file_pre, file_ext] = fileparts(input_pathfile);
file_pattern = [file_pre, '-*', file_ext]; % include hyphen, which is more specific to single trace ABR CSV files.
selpath = fileparts(selpath_analyze); % get parent folder containing "analyze" folder

% get listing of all single trace CSV files with the file name prefix
parent_folder = fullfile(selpath, file_pattern);
dir_files = dir(parent_folder);
n_files = numel(dir_files);

% Throw warning message and abort if parent directory contains no single trace ABR CSV Files.
if isempty(dir_files) 
    disp(['Warning! No single trace ABR CSV files in parent folder: ', parent_folder])
    disp('-> Aborting single trace analysis of ABR in get_roc_innerprod_ABR.m, line 98.')
    disp('')
    % Assign null values to outputs and return
    auc_cache = [];
    ip_mean_cache = [];
    ip_ste_cache = []; 
    thresh_cache = [];
    d_prime_cache_rms = [];
    d_prime_cache_z = [];
    return
end

% Calculate inner product distributions (i.e., histograms)
ip_dists_fulltrace = cell(n_files, 1);
ip_dists_window = cell(n_files, 1);
amps = zeros(n_files, 1);
freqs = zeros(n_files, 1);
sel_timewindow_range = cell(n_files, 1);
ave_X_csv = cell(n_files, 1);
dist_wave1amp = zeros(n_files, N_SINGLE_TRACES/2); % matrix of single trial wave 1 amplitudes. 
dist_ampnoise = zeros(n_files, N_SINGLE_TRACES/2); % matrix of single trial wave 1 amplitudes. 
% dc_single_trace = zeros(n_files, N_SINGLE_TRACES/2);
for i = 1:n_files
    
    % Load next single trial ABR data (CSV file)
    this_filename = dir_files(i).name;
    disp(['Working on file ', num2str(i), ' out of ', num2str(n_files), ': ', this_filename])
    [X_csv, A_csv, freq_csv] = import_ABRcsv(this_filename, selpath);
    
    % merge single trace pairs with alternating polarity to cancel cochlear
    % microphonic heterogeneity
    X_csv = merge_singletraceABR_polarities(X_csv); % 1026 -> 513 single traces. Matrix of size (SAMPLES, m_traces/2)
    
%     % 12-10-21 calculate DC shift of single traces
%     dc_single_trace(i, :) = mean(X_csv, 1);
    
    % Optional: restrict inner product analysis to a specific time window
    if TRIM_ABR 
        % Get time window for specific frequency and stimulus level of
        % single trace ABR CSV file that is loaded
        all_stimulus_levels = wave1_latency(:, 4);
        is_level = all_stimulus_levels==A_csv;
        
        all_freq = wave1_latency(:, 3);
        is_freq = all_freq == freq_csv;
        is_freq_level = is_level & is_freq;
        
        % convert units of time window from time (ms) to samples
        this_wave1_latency = int64(wave1_latency(is_freq_level, 1:2)*SAMPLING_RATE); % 2 element vector, [peak latency, trough latency], in units of sample index
        
        % Calculate time window (or time points) to "window"
        % single trace ABR using the wave 1 peak and trough latencies
        % 12-8-21: use single points
%         sel_timewindow_range{i} = this_wave1_latency; % use one point for peak, and 1 point for trough
%         sel_timewindow_range = max([sel_timewindow(1), 1]):sel_timewindow(2); % ensure time window starts no earlier than 0 ms
        
        % Use one time window between peak and trough, with buffer space
        % before peak and after trough
        timewindow_start = max([wave1_latency(is_freq_level, 1) - TIME_BEFOREPEAK, 0]); % ensure time window starts no earlier than 0 ms
%         timewindow_end = timewindow_start + TIME_WINDOW_LENGTH; % set time window duration from peak latency, regardless of trough latency
        timewindow_end = wave1_latency(is_freq_level, 2) + TIME_BEFOREPEAK; 
        %convert units ms -> sample index
        timewindow_start = int64(timewindow_start*SAMPLING_RATE);
        timewindow_end = int64(timewindow_end*SAMPLING_RATE);
        if timewindow_start==0
            timewindow_start = 1;
        end
        if timewindow_end==0
            timewindow_end = max([2, timewindow_start]);
        end
        sel_timewindow_range{i} = timewindow_start:timewindow_end; % use time window between peak and trough latency

%         % 12-22-21: Use two time windows, one centered around peak, second
%         % centered around trough
%         timewindow1_start = max([wave1_latency(is_freq_level, 1)-TIME_AROUNDPEAK/2, 0]); % ensure time window starts no earlier than 0 ms
%         timewindow1_end = timewindow1_start + TIME_AROUNDPEAK;
%         timewindow2_start = wave1_latency(is_freq_level, 2)-TIME_AROUNDPEAK/2;
%         timewindow2_end = timewindow2_start + TIME_AROUNDPEAK;
%         %convert units ms -> sample index
%         timewindow1_start = int64(timewindow1_start*SAMPLING_RATE);
%         timewindow1_end = int64(timewindow1_end*SAMPLING_RATE);
%         timewindow2_start = int64(timewindow2_start*SAMPLING_RATE);
%         timewindow2_end = int64(timewindow2_end*SAMPLING_RATE);
%         sel_timewindow_range{i} = [timewindow1_start:timewindow1_end, timewindow2_start:timewindow2_end]; % use 2 time windows around peak and trough latencies
        X_csv_window = X_csv(sel_timewindow_range{i}, :); % use time window(s) based on peak and trough latencies
%         X_csv_window = X_csv(timewindow1_start:timewindow2_start, :); % use time window between peak and trough 
        
        % Modified inner product method: 
        % Removes diagonal of matrix of inner products of single traces with
        % each other (N-1), before calculating mean inner product for single trace.
%         ip_dists_window{i} = get_innerprod_matrix_ABR(X_csv_window);
        
        % 12-29-21: commented original inner product distribution
        % calculation in line below.
        ip_dists_window{i} = get_innerprod_matrix_ABR(X_csv_window)./size(X_csv_window, 1); % divide inner product by number of samples in each product
        % 12-29-21: calculate wave 1 amp using a few points around peak and
        % trough
%         timewindow1_range = timewindow1_start:timewindow1_end;
%         timewindow2_range = timewindow2_start:timewindow2_end;
%         ip_dists_window{i} = get_wave1amp_singletrace(X_csv, timewindow1_range, timewindow2_range);
        % 12-31-21: calculate noise from each single trace, using random
        % point at start of each trace
        num_samples = size(X_csv, 1);
        X = (1/SAMPLING_RATE) * (1:num_samples); % time in ms
        this_distnoise = X_csv(1, :) - interp1(X, X_csv, 0.45); % noise in single trace
        dist_ampnoise(i, :) = reshape(this_distnoise, [length(this_distnoise), 1]); % make column vector
%         ip_dists_window{i} = X_csv(1, :) - interp1(X, X_csv, 0.45); % noise in single trace
%         ip_dists_window{i} = reshape(ip_dists_window{i}, [length(ip_dists_window{i}), 1]); % make column vector

%         ip_dists_window{i} = mean(X_csv_window(1, :), 2) - mean(X_csv_window(2, :), 2); % 12/9/21 - calculate wave 1 amp
    
        % 12-15-21: Calculate wave 1 amplitude of each single trace, to
        % calculate AUC and significance of difference between pre and post
        % 1 week level functions for wave 1 amplitude.
        num_samples = size(X_csv, 1);
        X = (1/SAMPLING_RATE) * (1:num_samples); % time in ms
        dist_wave1amp(i, :) = interp1(X, X_csv, wave1_latency(is_freq_level, 1)) - interp1(X, X_csv, wave1_latency(is_freq_level, 2)); % wave 1 amp = peak minus trough
    
    end
    
%     % Inner product method
%     ip_dists{i} = get_innerprod_raw_ABR(X_csv);
    
    % Modified inner product method: 
    % Removes diagonal of matrix of inner products of single traces with
    % each other (N-1), before calculating mean inner product for single trace.
    ip_dists_fulltrace{i} = get_innerprod_matrix_ABR(X_csv);
    
    amps(i) = A_csv;
    freqs(i) = freq_csv;
    
    % 12/9/21 - plot average ABRs, calculated from single traces, for
    % debugging why single trace window ABR is giving odd results (can't
    % reproduce wave 1 amp calculations)
    ave_X_csv{i} = mean(X_csv, 2);

end


%% Plot and analyze signal detection sensitivity from inner product distributions
freqs_unique = unique(freqs);
n_freq = length(freqs_unique);

amps_unique = unique(amps);
num_levels = length(amps_unique);

% Iterate over full trace analysis and windowed trace analysis, if
% calculated
if TRIM_ABR 
    all_ip_dists = {ip_dists_fulltrace, ip_dists_window};
else
    all_ip_dists = {ip_dists_fulltrace};
end
n_windows = numel(all_ip_dists);
signal_metrics_thresholds = cell(n_windows, 1); % output cache of signal detection sensitivity metrics and thresholds

if PLOT_FIGURES_DEFAULT
    % Obtain y axis scale for plotting average ABR
    ylim_max = [-1200, 1200];
    % Check to make sure bounds of plot are within y limits
    ymax_data = max(cellfun(@max, ave_X_csv));
    ymin_data = min(cellfun(@min, ave_X_csv));
    if ymax_data > ylim_max(2)
        ylim_max(2) = ymax_data;
    end
    if ymin_data < ylim_max(1)
        ylim_max(1) = ymin_data;
    end
end

for w = 1:n_windows
    % Load inner product distributions for full or windowed trace
    ip_dists = all_ip_dists{w};

    if PLOT_FIGURES_DEFAULT
        % Create tiled chart for displaying histograms
        fig = figure;
        t = tiledlayout(num_levels, n_freq, 'TileIndexing', 'columnmajor');         
        t.TileSpacing = 'tight';
        t.Padding = 'tight';
        
%         % 12/9/21 Create tiled chart for average ABRs, calculated from
%         % single traces, to validate single trial code
%         fig = figure;
%         t2 = tiledlayout(num_levels, n_freq, 'TileIndexing', 'columnmajor');         
%         t2.TileSpacing = 'tight';
%         t2.Padding = 'tight';
    end
    auc_cache = zeros(n_freq, num_levels);
    auc_cache_wave1 = zeros(n_freq, num_levels);
    ip_mean_cache = zeros(n_freq, num_levels);
    ip_ste_cache = zeros(n_freq, num_levels); % standard error of IP mean
    ip_std_cache = zeros(n_freq, num_levels); % standard deviation of IP distribution
    pval_cache = zeros(n_freq, num_levels);
    thresh_cache = ones(n_freq, 1)*MAX_THRESHOLD; % estimated ABR threshold per frequency
    thresh_criterion = ones(n_freq, 1)*MAX_THRESHOLD; % estimated ABR threshold per frequency
    cutoff = zeros(n_freq, 1);
    d_prime_cache_rms = zeros(n_freq, num_levels);
    thresh_d_prime_rms = ones(n_freq, 1)*MAX_THRESHOLD;
    d_prime_cache_z = zeros(n_freq, num_levels);
    thresh_d_prime_z = ones(n_freq, 1)*MAX_THRESHOLD;
    dist_wave1amp_cache = cell(n_freq, num_levels);
    dist_innerprod_cache = cell(n_freq, num_levels);

    % Calculate AUC, d', and mean values for each inner product histogram, 
    % iterating over frequencies and stimulus levels). These are the 
    % candidate metrics of signal detection sensitivity.
    for j = 1:n_freq
        % Iterate through frequencies
        this_freq = freqs_unique(j);

        % get data for this frequency only
        is_freq = freqs==this_freq;
        these_amps = amps(is_freq);
        these_ips = ip_dists(is_freq);
        these_ave_X_csv = ave_X_csv(is_freq);

        % Sort data in decreasing amplitude order
        [A_desc, ind_A_desc] = sort(these_amps, 'descend');
        ip_desc = these_ips(ind_A_desc);
        these_ave_X_csv_desc = these_ave_X_csv(ind_A_desc);
        dB = cell(num_levels, 1);
        dB(:) = {'dB (nV)'};
        newYlabels=cellfun(@(x,y) [x  ' ' y], cellstr(num2str(A_desc)), dB, 'un', 0);
 
        % Sort time points of wave 1 peak and trough if input
        if TRIM_ABR
            these_timewindow = sel_timewindow_range(is_freq);
            these_timewindow_desc = these_timewindow(ind_A_desc);
        end
        
        % initialize variables for skipping missing single trace ABR files
        missing_level_ind = NaN; % index of stimulus level if missing from set of single trace files at one frequency
        k_missing = 0; % counter for number of single trace ABR files are missing at frequency

        % check that single trace ABR files are available at all levels at this
        % frequency. If a file at one stimulus level is missing, label as NaN.
        this_amps_unique = unique(these_amps);
        this_num_levels = length(this_amps_unique);
        if this_num_levels ~= num_levels
            % Identify missing single trace ABR level(s) at this frequency
            [A_desc_all, ~] = sort(amps_unique, 'descend');
            [missing_levels, missing_level_ind] = setdiff(A_desc_all, A_desc, 'stable');
            disp(['Warning: missing single trace ABR file for stimulus level ', num2str(missing_levels), ' dB'])
        end

        % Added 11-12-21: Calculate ABR threshold using criterion cutoff value
        % set at noise inner product mean plus 3 standard deviations (check if
        % 3 is right constant)
        CUTOFF_MULT_STD = 2;
        cutoff(j) = mean(ip_desc{end}) + CUTOFF_MULT_STD*std(ip_desc{end});

        % Iterate through stimulus levels in descending order, skipping missing
        % levels
        is_thresh_cache = zeros(this_num_levels, 1);
        for k_all = 1:num_levels
            % Fill NaN values for signal detection metrics if single traces were not measured at
            % this stimulus level
            if any(missing_level_ind==k_all) 
                pval_cache(j, k_all) = NaN;
                auc_cache(j, k_all) = NaN;
                auc_cache_wave1(j, k_all) = NaN;
                ip_mean_cache(j, k_all) = NaN;
                d_prime_cache_rms(j, k_all) = NaN;
                d_prime_cache_z(j, k_all) = NaN;
                dist_wave1amp_cache{j, k_all} = NaN;
                dist_innerprod_cache{j, k_all} = NaN;
                k_missing = k_missing + 1;
                continue
            end
            k = k_all - k_missing; % drop index of stimulus level back if single trace ABR for a stimulus level is missing.

            % Fill NaN values for signal detection metrics if inner product distribution contains NaN values
            if any(isnan(ip_desc{k})) || any(isnan(ip_desc{end})) 
                pval_cache(j, k_all) = NaN;
                auc_cache(j, k_all) = NaN;
                auc_cache_wave1(j, k_all) = NaN;
                ip_mean_cache(j, k_all) = NaN;
                d_prime_cache_rms(j, k_all) = NaN;
                d_prime_cache_z(j, k_all) = NaN;
                dist_wave1amp_cache{j, k_all} = NaN;
                dist_innerprod_cache{j, k_all} = NaN;
                k_missing = k_missing + 1;
                continue
            end
            
            % calculate p-val
            p = ranksum(ip_desc{k}, ip_desc{end}, 'tail', 'right');
            pval_cache(j, k_all) = p;
            is_thresh = p<0.05;
            is_thresh_cache(k) = is_thresh;

            % Plot histograms
            if PLOT_FIGURES_DEFAULT
                % Plot histograms
                nexttile(t)
                histogram(ip_desc{k})
        %         xlim([-1, 1]) % if inner products are normalized 
        %         xlim([-2, 2]*10^-4) % if inner products are partially normalized 
%                 xlim([-0.5, 0.5]*10^-8) % if inner products are not normalized 
                if TRIM_ABR
%                     xlim([-5000, 10000]) % wave 1 amp dist limits
%                     frac_timepoints_inwindow = size(X_csv_window, 1) / size(X_csv, 1); 
%                     xlim(X_LIM_WHOLETRACE*frac_timepoints_inwindow) % maximum innerproduct value is proportional to number of samples in trace
                    xlim([-500000, 10000000]) % normalized window inner product dist limits
                else
                    xlim(X_LIM_WHOLETRACE) % if inner products are not normalized 
                end

                % Label x axis if last iteration
                if k_all==num_levels && j==n_freq
                    xlabel(t, 'Inner product (nV^2)')
                    ylabel(t, 'Frequency')
                end

                % Bold tile title if is at or above threshold
                if is_thresh
                    title([num2str(this_freq), ' Hz, ', num2str(A_desc(k)), ' dB'], 'fontweight','bold')
                else
                    title([num2str(this_freq), ' Hz, ', num2str(A_desc(k)), ' dB'], 'fontweight','normal')
                end
                
%                 % 12/9/21: Plot average ABR calculated with single trace data, for debugging
%                 nexttile(t2)
%                 plot(X, these_ave_X_csv_desc{k}, 'LineWidth', 3)
%                 ylim(ylim_max)
% 
%                 % show ylabels for first column only
%                 if j==1
%                     ylabel(newYlabels{k}, 'FontSize', 24)
%                     % rotate y label to make horizontal
%                     hYLabel = get(gca,'YLabel');
%                     set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right', 'FontSize', 24)
%                 end
% 
%                 if TRIM_ABR
%                     % Mark wave 1 peak and following trough
%                     lat_peak = these_timewindow_desc{k}(1);
%                     lat_trough = these_timewindow_desc{k}(2);
%                     peak_pt = [double(lat_peak)/SAMPLING_RATE, these_ave_X_csv_desc{k}(lat_peak)];
%                     trough_pt = [double(lat_trough)/SAMPLING_RATE, these_ave_X_csv_desc{k}(lat_trough)];
% %                     [peak_pt, trough_pt, amp, lat_peak, lat_trough] = these_timewindow_desc
% %                     if k==1
% %                         % Get wave 1 peak and following trough at highest stimulus
% %                         % level
% %                         [peak_pt, trough_pt, amp, lat_peak, lat_trough] = get_wave1_averageABR(X, y);
% %                     else
% %                         % Get wave 1 peak and following trough at lower stimulus
% %                         % level. Ensure peak latency does not decrease.
% %                         [peak_pt, trough_pt, amp, lat_peak, lat_trough] = get_wave1_averageABR(X, y, peak_pt(1), trough_pt(1));
% %                     end
%                     hold on
%                     CIRCLE_SIZE = 54; % default 36
% %                     % 12-14-21: calculate wave 1 peak and trough using
% %                     % original input to function, for debugging
% %                     count = k + (j-1)*num_levels;
% %                     peak_pt(1) = wave1_latency(count, 1);
% %                     trough_pt(1) = wave1_latency(count, 2);
% %                     peak_pt(2) = interp1(X, these_ave_X_csv_desc{k}, peak_pt(1));
% %                     trough_pt(2) = interp1(X, these_ave_X_csv_desc{k}, trough_pt(1));
%                     scatter(peak_pt(1), peak_pt(2), CIRCLE_SIZE, 'green', 'filled') % peak
%                     scatter(trough_pt(1), trough_pt(2), CIRCLE_SIZE, 'red', 'filled') % trough
%                     hold off
%     
%                     % Display wave 1 measurements in corner of plot
%                     xL=xlim;
%                     yL=ylim;
%                     amp = peak_pt(2) - trough_pt(2);
%                     str = sprintf('P-P_I = %.0f nV', amp);
% %                     if A_descending(i) > thresh_cache_liberman
% %                         text(0.99*xL(2), 0.99*yL(1), str,'HorizontalAlignment','right','VerticalAlignment','bottom', 'FontWeight', 'bold')
% % %                         is_abovenoise(i,1)=1;
% %                     else
%                     text(0.99*xL(2), 0.99*yL(1), str,'HorizontalAlignment','right','VerticalAlignment','bottom')
% %                         is_abovenoise(i,1)=0;
% %                     end
%                 end
% 
%                 % Remove extraneous axis tick marks and x-axis from all but bottom
%                 % tile
%                 set(gca,'box','off')
%                 if k~=num_levels
%                     set(gca,'xtick',[])
%                     set(gca,'xticklabel',[])
%                     set(gca,'XColor','none')
%                 end
% 
%                 % Remove y-axis labels from all but first column
%                 if j~=1
%                     set(gca,'yticklabel',[])
%                 end
% 
%                 % Set text size of current axis
%                 set(gca,'FontSize',24)
% 
%                 % Show frequency as title above top-most tile of each column
%                 if k==1
%                     title([num2str(this_freq), ' Hz'])
%                 end
%                 
%                 if j==n_freq && k_all==num_levels
%                     xlabel(t2, 'Time (ms)', 'FontSize', 24)
%                     t2.TileSpacing = 'none';
%                     t2.Padding = 'tight';
% 
%                     % Maximize figure window size
%                     fig.WindowState = 'maximized';
%                 end
            end

            % calculate AUC value at this level
            scores = vertcat(ip_desc{k}, ip_desc{end});
            num_pts = numel(scores);
            labels = zeros(num_pts, 1);
            labels(1:num_pts/2, 1) = 1;
            posclass = 1;
            
            % Fill NaN values for signal detection metrics if unable to
            % calculate AUC value due to NaN values in scores
            if any(isnan(scores))
                pval_cache(j, k_all) = NaN;
                auc_cache(j, k_all) = NaN;
                auc_cache_wave1(j, k_all) = NaN;
                ip_mean_cache(j, k_all) = NaN;
                d_prime_cache_rms(j, k_all) = NaN;
                d_prime_cache_z(j, k_all) = NaN;
                dist_wave1amp_cache{j, k_all} = NaN;
                dist_innerprod_cache{j, k_all} = NaN;
                continue
            end
            % calculate AUC without bootstrapping
            [~,~,~,AUC] = perfcurve(labels,scores,posclass);
            auc_cache(j, k_all) = AUC;

    %         %Bootstrapping assumes feature cutoff is fixed at each point in ROC, ie Tvals choice
    %         nboot = 5; %nboot = 500 => convergence to less than 10% error
    %         [~,~,~,AUC] = perfcurve(labels,scores,posclass,'NBoot',nboot,'BootType','stud','Tvals','All','Options',statset('UseParallel',false));
    %         auc_cache(j, k_all) = AUC(1);
    %         auc_lower_cache(j, k_all) = AUC(2);
    %         auc_upper_cache(j, k_all) = AUC(3);

            % Calculate mean inner product value of single traces with coherent
            % average
            ip_mean_cache(j, k_all) = mean(ip_desc{k});
            ip_std_cache(j, k_all) = std(ip_desc{k});
            ip_ste_cache(j, k_all) = std(ip_desc{k})/sqrt(length(ip_desc{k}));

            % Calculate d' signal detection sensitivity. 
            d_prime_cache_rms(j, k_all) = (mean(ip_desc{k}) - mean(ip_desc{end}))/rms([std(ip_desc{k}), std(ip_desc{end})]); % d' = u_signal - u_noise / rms(std_signal, std_noise)
            d_prime_cache_z(j, k_all) = (cutoff(j) - mean(ip_desc{end}))/std(ip_desc{end}) - (cutoff(j) - ip_mean_cache(j, k_all))/ip_std_cache(j, k_all); % d' = z_noise(criterion) - z_signal(criterion)
        
            % Cache distribution of wave 1 amplitude and inner product 
            % values
            this_dist = dist_wave1amp(is_freq, :); 
            this_ampnoise_dist = dist_ampnoise(is_freq, :); % 12-31-21: try calculating AUC using noise dist of points in same stimulus level single trace, at noisy points at start of each trace
            dist_wave1amp_cache{j, k_all} = this_dist(ind_A_desc(k), :);
            dist_innerprod_cache{j, k_all} = ip_desc{k};
            
            % Calculate AUC for wave 1 amplitude distribution across single
            % traces
            scores = horzcat(this_dist(ind_A_desc(k), :), this_dist(ind_A_desc(end), :));
%             scores = horzcat(this_dist(ind_A_desc(k), :),
%             this_ampnoise_dist(ind_A_desc(k), :)); % 12-31-21: use first point of
%             each trace to determine noise distribution
            num_pts = numel(scores);
            labels = zeros(num_pts, 1);
            labels(1:num_pts/2, 1) = 1;
            posclass = 1;
            % Fill NaN values for signal detection metrics if unable to
            % calculate AUC value due to NaN values in scores
            if any(isnan(scores))
                pval_cache(j, k_all) = NaN;
                auc_cache(j, k_all) = NaN;
                auc_cache_wave1(j, k_all) = NaN;
                ip_mean_cache(j, k_all) = NaN;
                d_prime_cache_rms(j, k_all) = NaN;
                d_prime_cache_z(j, k_all) = NaN;
                dist_wave1amp_cache{j, k_all} = NaN;
                dist_innerprod_cache{j, k_all} = NaN;
                continue
            end
            % calculate AUC without bootstrapping
            [~,~,~,AUC] = perfcurve(labels,scores,posclass);
            auc_cache_wave1(j, k_all) = AUC;
        end

        % Calculate threshold at this frequency using Wilcoxin rank-sum p-value
        ind_first_belowthresh = find(~is_thresh_cache, 1);
        ind_thresh = ind_first_belowthresh - 1;
        if ind_thresh >= 1
            thresh_cache(j) = A_desc(ind_thresh);
        end

        % Added 11-12-21: Calculate ABR threshold using criterion cutoff value
        % set at noise inner product mean plus 3 standard deviations (check if
        % 3 is right constant)
    %     CUTOFF_MULT_STD = 2;
    %     cutoff(j) = ip_mean_cache(j, end) + CUTOFF_MULT_STD*ip_std_cache(j, end);
        ind_thresh_cutoff = find(ip_mean_cache(j, :) > cutoff, 1, 'last');
        if ind_thresh_cutoff >= 1
    %         thresh_criterion(j) = A_desc(ind_thresh_cutoff);
            thresh_criterion(j) = xgiveny_avetemplate(cutoff(j), flip(ip_mean_cache(j, :)), flip(A_desc));
        end

        % Calculate threshold using d'
        CUTOFF_D_PRIME_RMS = 0.5;
        thresh_d_prime_rms(j) = xgiveny_avetemplate(CUTOFF_D_PRIME_RMS, flip(d_prime_cache_rms(j, :)), flip(A_desc));
        CUTOFF_D_PRIME_Z = 1;
        thresh_d_prime_z(j) = xgiveny_avetemplate(CUTOFF_D_PRIME_Z, flip(d_prime_cache_z(j, :)), flip(A_desc));
    
    end

    % Reshape cache of wave 1 amplitude and inner product values to
    % cell of size (nfreq, 1) to fit in one column
    dist_wave1amp_cache_temp = cell(n_freq, 1);
    dist_innerprod_cache_temp = cell(n_freq, 1);
    for zz=1:n_freq
        dist_wave1amp_cache_temp{zz} = dist_wave1amp_cache(zz, :);
        dist_innerprod_cache_temp{zz} = dist_innerprod_cache(zz, :);
    end
    dist_wave1amp_cache = dist_wave1amp_cache_temp;
    dist_innerprod_cache = dist_innerprod_cache_temp;
    
    % Cache signal detection sensitivity measures for this level of
    % windowing.
    signal_metrics_thresholds{w} = {auc_cache, ip_mean_cache, ip_ste_cache, thresh_cache, thresh_criterion, d_prime_cache_rms, thresh_d_prime_rms, d_prime_cache_z, thresh_d_prime_z, dist_wave1amp_cache, dist_innerprod_cache, auc_cache_wave1};


%     %% Plot stimulus vs AUC curves for each frequency
%     if ~PLOT_FIGURES_DEFAULT
%         continue
%     end
% 
%     figure
%     hold on
%     for kk = 1:n_freq
%         plot(A_desc, auc_cache(kk, :))
%     end
%     ylabel('AUC')
%     xlabel('dB')
%     legend(num2cell(string(freqs_unique)))
% 
%     %% Plot stimulus vs IP mean for each frequency
%     figure
%     hold on
%     for kkk = 1:n_freq
%         plot(A_desc, ip_mean_cache(kkk, :))
%         disp(['Criterion threshold @ ', num2str(freqs_unique(kkk)), ' Hz: ', num2str(thresh_criterion(kkk)), ' dB'])
%     end
%     ylabel('Inner product, mean (nV^2)')
%     xlabel('dB')
%     legend(num2cell(string(freqs_unique)))
% 
%     % Plot criterion 
%     yline(cutoff)

    % % Draw vertical line at threshold determine by rank sum of innerprod
    % % distribution
    % yL = get(gca,'YLim');
    % line([4 4],yL,'Color','r');
    % line([7 7],yL,'Color','r');

%     %% Plot stimulus vs d' (calculated using d' = u_signal - u_noise / rms(std_signal, std_noise)) for each frequency
%     figure
%     hold on
%     for kkk = 1:n_freq
%         plot(A_desc, d_prime_cache_rms(kkk, :))
%         disp(['D prime rms threshold @ ', num2str(freqs_unique(kkk)), ' Hz: ', num2str(thresh_d_prime_rms(kkk)), ' dB'])
%     end
%     ylabel('D prime (a.u.)')
%     xlabel('dB')
%     legend(num2cell(string(freqs_unique)))
% 
%     %% Plot stimulus vs d' (calculated using d' = z_noise(criterion) - z_signal(criterion)) for each frequency
%     figure
%     hold on
%     for kkk = 1:n_freq
%         plot(A_desc, d_prime_cache_z(kkk, :))
%         disp(['D prime z threshold @ ', num2str(freqs_unique(kkk)), ' Hz: ', num2str(thresh_d_prime_z(kkk)), ' dB'])
%     end
%     ylabel('D prime (a.u.)')
%     xlabel('dB')
%     legend(num2cell(string(freqs_unique)))
end

end