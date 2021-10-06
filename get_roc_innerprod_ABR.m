function [auc_cache, ip_mean_cache, ip_ste_cache, thresh_cache, auc_lower_cache, auc_upper_cache] = get_roc_innerprod_ABR(varargin)
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
%       - varargin{3}: optional specification of start and end time points
%       to perform single trace analysis on. 2 element vector [start, end]
%       times in milliseconds.
%
% Outputs:
%       - auc_cache:    ROC-AUC values for each set of single ABR traces.
%                       Matrix of size (n_freq, num_levels).
%       - ip_mean_cache: Mean inner product value of all single traces'
%                       inner products with coherent average at same stimulus 
%                       level. Matrix of size (n_freq, num_levels). 
%       - ip_ste_cache: Standard error of inner product mean. Matrix of size (n_freq, num_levels). 
%       - thresh_cache: Estimated ABR thresholds, based on Wilcoxin
%                       sum-rank test p-val < 0.5 for inner product
%                       distribution with noise distribution. Vector of length (n_freq). In order of freqs_unique. 
%       - auc_lower_cache:    ROC-AUC lower bound confidence estimate based on bootstrapping.
%                       Matrix of size (n_freq, num_levels).
%       - auc_upper_cache:    ROC-AUC upper bound confidence estimate based on bootstrapping.
%                       Matrix of size (n_freq, num_levels).
%
% Example usage: 
% [auc_cache, ip_mean_cache, thresh_cache, A_desc, freqs_unique] = get_roc_innerprod_ABR;
%
% George Liu
% Last edited 9/8/21 
% - added 3rd input option to allow innerprod analysis on selected time window of ABR trace.
% - added criteria that ABR threshold require p-val < 0.5 at AND for all
% stimulus levels above threshold.
% - added bootstrapping to estimate AUC lower and upper confidence intervals for standard errors
% - calculate standard error for inner product mean
% 
% 9/1/21 - save NaN values for missing single trace ABR files
% (eg if first set of single trace files were accidentally during saving of average ABR files).
%
% Dependencies: get_innerprod_raw_ABR.m, get_innerprod_matrix_ABR.m, merge_singletraceABR_polarities.m

% [filename, path] = uigetfile('*.csv');
% [X_csv, A_csv, freq_csv] = import_ABRcsv(filename, path);
% this_ip_lone = get_innerprod_raw_ABR(X_csv);
% figure, histogram(this_ip_lone)

%% Constants
MAX_THRESHOLD = 95; % dB; arbitrary max threshold if ABRs are noisy at all stimuli levels
SAMPLING_RATE = 200; % samples / msec. Of note, average CSV files have different sampling rate b/c they compress number of samples 1952 -> 243.
TRIM_ABR = false;

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
    trim_abr_times = varargin{3};
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

% Throw error if parent directory contains no single trace ABR CSV Files.
if isempty(dir_files) 
    error(['Error ocurred. No single trace ABR CSV files in folder: ', parent_folder])
end

ip_dists = cell(n_files, 1);
amps = zeros(n_files, 1);
freqs = zeros(n_files, 1);
for i = 1:n_files
    this_filename = dir_files(i).name;
    disp(['Working on file ', num2str(i), ' out of ', num2str(n_files), ': ', this_filename])
    [X_csv, A_csv, freq_csv] = import_ABRcsv(this_filename, selpath);
    
    % merge single trace pairs with alternating polarity to cancel cochlear
    % microphonic heterogeneity
    X_csv = merge_singletraceABR_polarities(X_csv); % 1026 -> 513 single traces
    
    % Optional: restrict inner product analysis to a specific time window
    if TRIM_ABR 
        sel_timewindow = int64(trim_abr_times*SAMPLING_RATE);
        sel_timewindow_range = sel_timewindow(1):sel_timewindow(2);
        X_csv = X_csv(sel_timewindow_range, :);
    end
    
%     % Inner product method
%     ip_dists{i} = get_innerprod_raw_ABR(X_csv);
    
    % Modified inner product method: 
    % Removes diagonal of matrix of inner products of single traces with
    % each other (N-1), before calculating mean inner product for single trace.
    ip_dists{i} = get_innerprod_matrix_ABR(X_csv);
    
    amps(i) = A_csv;
    freqs(i) = freq_csv;
end

%%
freqs_unique = unique(freqs);
n_freq = length(freqs_unique);

amps_unique = unique(amps);
num_levels = length(amps_unique);

if PLOT_FIGURES_DEFAULT
    figure
    t = tiledlayout(num_levels, n_freq, 'TileIndexing', 'columnmajor');
end
auc_cache = zeros(n_freq, num_levels);
auc_lower_cache = zeros(n_freq, num_levels); % lower bound estimate for AUC standard error
auc_upper_cache = zeros(n_freq, num_levels);
ip_mean_cache = zeros(n_freq, num_levels);
ip_ste_cache = zeros(n_freq, num_levels); % standard error of IP mean
pval_cache = zeros(n_freq, num_levels);
thresh_cache = ones(n_freq, 1)*MAX_THRESHOLD; % estimated ABR threshold per frequency

for j = 1:n_freq
    this_freq = freqs_unique(j);
    
    % get data for this frequency only
    is_freq = freqs==this_freq;
    these_amps = amps(is_freq);
    these_ips = ip_dists(is_freq);
    
    % Sort data in decreasing amplitude order
    [A_desc, ind_A_desc] = sort(these_amps, 'descend');
    ip_desc = these_ips(ind_A_desc);
    
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
    
    % Iterate through stimulus levels in descending order, skipping missing
    % levels
    is_thresh_cache = zeros(this_num_levels, 1);
    for k_all = 1:num_levels
        % Fill NaN values if this level is missing
        if any(missing_level_ind==k_all)
            pval_cache(j, k_all) = NaN;
            auc_cache(j, k_all) = NaN;
            ip_mean_cache(j, k_all) = NaN;
            k_missing = k_missing + 1;
            continue
        end
        k = k_all - k_missing; % drop index of stimulus level back if single trace ABR for a stimulus level is missing.
        
        % calculate p-val
        p = ranksum(ip_desc{k}, ip_desc{end}, 'tail', 'right');
        pval_cache(j, k_all) = p;
        is_thresh = p<0.05;
        is_thresh_cache(k) = is_thresh;
        
        % Plot histograms
        if PLOT_FIGURES_DEFAULT
            nexttile
            histogram(ip_desc{k})
    %         xlim([-1, 1]) % if inner products are normalized 
    %         xlim([-2, 2]*10^-4) % if inner products are partially normalized 
            xlim([-0.5, 0.5]*10^-8) % if inner products are not normalized 
            
            % Label x axis if last iteration
            if k_all==num_levels && j==n_freq
                xlabel(t, 'Inner product (V^2)')
                ylabel(t, 'Frequency')
            end
            
            % Bold tile title if is at or above threshold
            if is_thresh
                title([num2str(this_freq), ' Hz, ', num2str(A_desc(k)), ' dB'], 'fontweight','bold')
            else
                title([num2str(this_freq), ' Hz, ', num2str(A_desc(k)), ' dB'], 'fontweight','normal')
            end
        end
        
        % calculate AUC value at this level
        scores = vertcat(ip_desc{k}, ip_desc{end});
        num_pts = numel(scores);
        labels = zeros(num_pts, 1);
        labels(1:num_pts/2, 1) = 1;
        posclass = 1;
%         % calculate AUC without bootstrapping
%         [~,~,~,AUC] = perfcurve(labels,scores,posclass);
%         auc_cache(j, k_all) = AUC;
        
        %Bootstrapping assumes feature cutoff is fixed at each point in ROC, ie Tvals choice
        nboot = 5; %nboot = 500 => convergence to less than 10% error
        [~,~,~,AUC] = perfcurve(labels,scores,posclass,'NBoot',nboot,'BootType','stud','Tvals','All','Options',statset('UseParallel',false));
        auc_cache(j, k_all) = AUC(1);
        auc_lower_cache(j, k_all) = AUC(2);
        auc_upper_cache(j, k_all) = AUC(3);
        
        % Calculate mean inner product value of single traces with coherent
        % average
        ip_mean_cache(j, k_all) = mean(ip_desc{k});
        ip_ste_cache(j, k_all) = std(ip_desc{k})/sqrt(length(ip_desc{k}));
        
    end
    
    % Calculate threshold at this frequency
    ind_first_belowthresh = find(~is_thresh_cache, 1);
    ind_thresh = ind_first_belowthresh - 1;
    if ind_thresh>=1
        thresh_cache(j) = A_desc(ind_thresh);
    end
    
end


%% Plot stimulus vs AUC curves for each frequency
if ~PLOT_FIGURES_DEFAULT
    return
end

figure
hold on
for kk = 1:n_freq
    plot(A_desc_all, auc_cache(kk, :))
end
ylabel('AUC')
xlabel('dB')
legend(num2cell(string(freqs_unique)))

%% Plot stimulus vs IP mean for each frequency
figure
hold on
for kkk = 1:n_freq
    plot(A_desc_all, ip_mean_cache(kkk, :))
end
ylabel('Inner product, mean (V^2)')
xlabel('dB')
legend(num2cell(string(freqs_unique)))

% % Draw vertical line at threshold determine by rank sum of innerprod
% % distribution
% yL = get(gca,'YLim');
% line([4 4],yL,'Color','r');
% line([7 7],yL,'Color','r');

end