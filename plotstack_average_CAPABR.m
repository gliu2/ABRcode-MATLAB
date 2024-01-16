function plotstack_average_CAPABR(varargin)
% Plot a stack of average CAPs and ABRs using single traces, one stack per frequency. Saves
% plots in multiple image formats.
%
% Default is for user to select file in dialog box when no input parameters
% are specified.
%
% Optional input syntax for specifying input file: plotstack_averageABR(path, filename)
%
% TODO: 
% [ ] Some TDT systems save average ABR CSV file in units of V instead of
% nV. Make import_averageABRcsv.m check units of ABR voltages and convert
% values to nV if necessary. (Added 10/6/21).
%
% 6/30/2023 George Liu
% Dependencies: get_wave1_averageCAP.m, get_wave1_averageABR.m, get_fft_CAPABR.m, 
%           notch_filter.m, filter_lowpass_blackman.m, filter_highpass_blackman.m, 
%           import_ABRcsv.m, merge_singletraceABR_polarities, savefig_multiformat.m

%% Constants
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures
% SAMPLE_PERIOD = 40.96; % sample period of 40.96 us/sample in average trace
SAMPLE_PERIOD = 5.12; % sample period of 5.12 us/sample in single trace
SAMPLE_PERIOD_MS = SAMPLE_PERIOD/1000; % sample period in ms

NUM_CHANNELS = 2; % 1 = ABR, 2 = CAP
CHANNEL_KEY = {'ABR', 'CAP'};
CHANNEL_YLIM = {[-1000, 1000], [-1000, 1000]};

NUM_LEVELS = 10;
% N_FREQ = 2; % 16 kHz, 32 kHz

N_ORDER = 700; % Half mainlobe width (-3 dB, approximation of cutoff frequency) of 226.484 Hz
% F_HIGHPASS = 100; % Hz
% F_LOWPASS = 3000; % Hz

% Default settings for when running function directly
FILTER_HIGHPASS_ON = 1; % Apply high pass filter with Blackman window of N_ORDER length. Use if no preacquisition high pass filtering was done.
FILTER_LOWPASS_ON = 1; % Apply low pass filter with cutoff of 3.052 kHz (Blackman window of order 100). Use if no preacquisition low pass filtering was done
FILTER_NOTCH_ON = 0;

PLOT_FIGURES = 1; % Hard code logical for whether or not to plot figures 
PLOT_WAVE1 = 1; % Plot overlay circles on peak and trough of wave 1
SAVE_FIGURES = 1; % Hard code logical for whether or not to save figures 

% Load average trace data
if nargin == 2
    path = varargin{1};
    filename = varargin{2};
elseif nargin == 0
    [filename, path] = uigetfile('*.csv');
elseif nargin == 1
    [path, name, ext] = fileparts(varargin{1});
    filename = [name, ext];
elseif nargin == 3
    path = varargin{1};
    filename = varargin{2};
    FILTER_HIGHPASS_ON = varargin{3};
elseif nargin == 4
    path = varargin{1};
    filename = varargin{2};
    FILTER_HIGHPASS_ON = varargin{3};
    PLOT_FIGURES = varargin{4};
elseif nargin == 5
    path = varargin{1};
    filename = varargin{2};
    FILTER_HIGHPASS_ON = varargin{3};
    PLOT_FIGURES = varargin{4};
    SAVE_FIGURES = varargin{5};
else
    disp('Warning: number of input arguments is not valid!')
    [filename, path] = uigetfile('*.csv');
end

% [filename,path] = uigetfile('*.csv'); % Uncomment if running as script
% rather than function
disp(['Opening ', fullfile(path, filename)])

% Analyze metadata in single trial ABR/CAP filename
% Filename example: 20230628_tmie2_ABRCAP_pre-0-55-2-1.csv
% group 0
% SGI 55
% channel 2
% run (?) 1
dash_loc = strfind(filename, '-');
dot_loc = strfind(filename, '.');

num_dashes = length(dash_loc);
group_ind = dash_loc(num_dashes - 3) + 1 : dash_loc(num_dashes - 2) - 1;
SGI_ind = dash_loc(num_dashes - 2) + 1 : dash_loc(num_dashes - 1) - 1;
channel_ind = dash_loc(num_dashes - 1) + 1 : dash_loc(num_dashes) - 1;
run_ind = dash_loc(num_dashes) + 1 : dot_loc - 1;

% Metadata of selected file
group = str2double(filename(group_ind));
% SGI = str2double(filename(SGI_ind));
% channel = str2double(filename(channel_ind)); % 1 = ABR, 2 = CAP
run = str2double(filename(run_ind));

% Obtain all single trial ABR/CAP files from same "run" / timepoint, which
% have the same filename stem.
filename_stem = filename(1 : dash_loc(end-3) - 1);
listing = dir(fullfile(path, [filename_stem, '*.csv']));

allnames = {listing.name}'; % column cell array of filenames

extract = @(C, k) cellfun(@(c) str2double((c(k))), C) ; % cell function to extract k'th element from each cell in cell array C.
all_groups = extract(allnames, group_ind);
all_SGI = extract(allnames, SGI_ind);
all_channels = extract(allnames, channel_ind); % 1 = ABR, 2 = CAP
all_runs = extract(allnames, run_ind);

for i = 1:NUM_CHANNELS
% for i = 2
    disp(['Working on channel ', num2str(i), ' out of ', num2str(NUM_CHANNELS), ': ', CHANNEL_KEY{i}])
    % Obtain ABR/CAP traces, ordered by SGI (same group and run as selected trace)
    is_data = all_channels == i;
    is_group = all_groups == group;
    is_run = all_runs == run;

    data_filenames = allnames(is_data & is_group & is_run);
    [~, sgi_order] = sort(all_SGI(is_data & is_group & is_run));
    data_filenames_ordered = data_filenames(sgi_order);
    
    num_data = length(data_filenames_ordered);
%     num_levels = num_data / N_FREQ; % assume only 2 tone pip frequencies, 16 and 32 kHz
    n_freq = num_data / NUM_LEVELS; % assume only 10 stimulus levels, 0, 10, ..., 90 dB
    
    if PLOT_FIGURES
        fig = figure;
        t = tiledlayout(NUM_LEVELS, n_freq, 'TileIndexing', 'columnmajor');
    end
        
    % Collect all average ABR traces from single trial ABR recordings
    y_avg_cache = cell(num_data, 1);
    y_avg_cm_cache = cell(num_data, 1); % microphonic
    A_csv_cache = zeros(num_data, 1);
    freq_csv_cache = zeros(num_data, 1);
    X_csv_cache = cell(num_data, 1);
    wave1peak_cache = cell(num_data, 1);
    wave1trough_cache = cell(num_data, 1);
    wave1amp_cache = zeros(num_data, 1);
    wave1lat_peak_cache = zeros(num_data, 1);
    wave1lat_trough_cache = zeros(num_data, 1);
    dist_wave1amp = cell(num_data, 1);
    wave1std = zeros(num_data, 1);

    for j = 1:num_data
        this_filename = data_filenames_ordered{j};
        disp(['Working on file ', num2str(j), ' out of ', num2str(num_data), ': ', this_filename])
        [X_csv, A_csv, freq_csv] = import_ABRcsv(this_filename, path);
        
        % Store cochlear microphonic prior to post-acquisition (low pass)
        % filtering
        y_cm = mean(merge_singletraceABR_polarities_CM(X_csv), 2)'; % Row vector
        
        % Apply post-acquisition filtering if pre-acquisition filtering
        % not done
        if FILTER_HIGHPASS_ON 
%             X_csv = filter_highpass_movingavg(X_csv, N_ORDER); % Moving average high pass filter
            X_csv = filter_highpass_blackman(X_csv, N_ORDER); % Blackman high pass filter
%              X_csv = filter_singletraceABR(X_csv, N_ORDER, F_HIGHPASS); %
%              Blackman high pass filter - stopped using because shifts
%              curve
%             X_csv = filter_singletraceABR(X_csv, N_ORDER, F_HIGHPASS,
%             F_LOWPASS); % Blackman bandpass filter - stopped using because shifts
%              curve
        end
        if FILTER_LOWPASS_ON 
            X_csv = filter_lowpass_blackman(X_csv, 150); % Blackman low pass filter with cutoff of 1.05 kHz (Blackman window of order 150)
            X_csv = filter_lowpass_blackman(X_csv, 150);  % Apply twice to suppress cycle (and alias) noise more
%             X_csv = filter_lowpass_blackman(X_csv, 86); % Blackman low pass filter with cutoff of 1812 Hz, for anti-aliasing and removal of occasional cycle noise at 1.9 kHz
        end
        if FILTER_NOTCH_ON 
            % Apply notch filter to remove 1901.1 Hz noise
            X_csv = notch_filter(X_csv);
        end
    
        % merge single trace pairs with alternating polarity to cancel cochlear
        % microphonic heterogeneity
        y = merge_singletraceABR_polarities(X_csv); % 1026 -> 513 single traces. Matrix of size (SAMPLES, m_traces/2)
        y_avg = mean(y, 2);
                   
        % Calculate number of samples, though this actually only needs to
        % be calculated once....
        if j==1
            num_samples = size(y, 1);
    %         n_traces = size(y, 2);
            X = SAMPLE_PERIOD_MS * (1:num_samples); % time in ms
        end      
        
        % Calculate wave 1 amplitude and latency
        if mod(j, NUM_LEVELS)==1
            % Get wave 1 peak and following trough at highest stimulus
            % level
            if i==1 % ABR
                [peak_pt, trough_pt, amp, lat_peak, lat_trough] = get_wave1_averageABR(X, y_avg);
            else % CAP
                [peak_pt, trough_pt, amp, lat_peak, lat_trough] = get_wave1_averageCAP(X, y_avg);
            end
        else
            % Get wave 1 peak and following trough at lower stimulus
            % level. Ensure peak latency does not decrease.
            if i==1 % ABR
                [peak_pt, trough_pt, amp, lat_peak, lat_trough] = get_wave1_averageABR(X, y_avg, peak_pt(1), trough_pt(1));
            else % CAP
                [peak_pt, trough_pt, amp, lat_peak, lat_trough] = get_wave1_averageCAP(X, y_avg, peak_pt(1), trough_pt(1));
            end
        end
        
        % Calculate wave 1 amplitude variability across single trials
        dist_wave1amp{j} = interp1(X, X_csv, lat_peak) - interp1(X, X_csv, lat_trough); % wave 1 amp = peak minus trough
        if i==2 % CAP
            dist_wave1amp{j} = -1*dist_wave1amp{j}; % if CAP, invert amplitude so it is generally positive
        end
        wave1std(j) = std(dist_wave1amp{j}); % standard deviation of wave 1 amplitude across single trials
              
        % Calculate CM for 90 dB, 8 kHz stimulus
        if j==1
            [f, P1, ~] = get_fft_CAPABR(y_cm);
            amp_cm_8khz90db = interp1(f, P1, 8000);
        end
        
        % Cache variables
        y_avg_cache{j, 1} = y_avg'; % make row vector
        y_avg_cm_cache{j, 1} = y_cm; % Row vector
        A_csv_cache(j, 1) = A_csv;
        freq_csv_cache(j, 1) = freq_csv;
        X_csv_cache{j, 1} = X_csv;
        wave1peak_cache{j} = peak_pt;
        wave1trough_cache{j} = trough_pt;
        wave1amp_cache(j) = amp;
        wave1lat_peak_cache(j) = lat_peak;
        wave1lat_trough_cache(j) = lat_trough;
        
    end
        
    % Make sure bounds of plot are within y limits
    ylim_max = CHANNEL_YLIM{i};
    ymax_data = max(cell2mat(y_avg_cache), [], "all");
    ymin_data = min(cell2mat(y_avg_cache), [], "all");
    if ymax_data > ylim_max(2)
        ylim_max(2) = ymax_data;
    end
    if ymin_data < ylim_max(1)
        ylim_max(1) = ymin_data;
    end
            
    if PLOT_FIGURES
        for j = 1:num_data
            y_avg = y_avg_cache{j};
            A_csv = A_csv_cache(j);
            freq_csv = freq_csv_cache(j);
            
            nexttile(t)
            plot(X, y_avg, 'LineWidth', 3)
            ylim(ylim_max)
            
            % Show CAP wave 1 peak and trough
            if PLOT_WAVE1 
                hold on
                CIRCLE_SIZE = 54; % default 36
                peak_pt = wave1peak_cache{j};
                trough_pt = wave1trough_cache{j};
                amp = wave1amp_cache(j);
                
                scatter(peak_pt(1), peak_pt(2), CIRCLE_SIZE, 'green', 'filled') % peak
                scatter(trough_pt(1), trough_pt(2), CIRCLE_SIZE, 'red', 'filled') % trough
%                 disp(['Wave 1 amp (nv): ', num2str(amp)])
                
                % Display wave 1 measurements in corner of plot
                xL=xlim;
                yL=ylim;
                str = sprintf('P-P_I = %.0f nV', amp);
                text(0.99*xL(2), 0.99*yL(1), str,'HorizontalAlignment','right','VerticalAlignment','bottom')
                hold off
            end
                        
            % show ylabels for first column only
            if j/NUM_LEVELS <= 1
                ylabel([num2str(A_csv), ' dB (nV)'], 'FontSize', 24)
                % rotate y label to make horizontal
                hYLabel = get(gca,'YLabel');
                set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right', 'FontSize', 24)
            end

            % Adjust plot appearance
            % Remove extraneous axis tick marks and x-axis from all but bottom
            % tile
            set(gca,'box','off')
            if mod(j, NUM_LEVELS) ~= 0 % 
                set(gca,'xtick',[])
                set(gca,'xticklabel',[])
                set(gca,'XColor','none')
            end

            % Remove y-axis labels from all but first column
            if j/NUM_LEVELS > 1
                set(gca,'yticklabel', [])
            end

            % Set text size of current axis
            set(gca,'FontSize',24)

            % Show frequency as title above top-most tile of each column
            if mod(j, NUM_LEVELS)==1
                title([num2str(freq_csv), ' Hz']) 
            end
        end
        
        % Format figure
        xlabel(t, 'Time (ms)', 'FontSize', 24)
    %     ylabel(t, 'Voltage (nV)', 'FontSize', 24)
        title(t, CHANNEL_KEY{i}, 'FontSize', 24)

        t.TileSpacing = 'none';
        t.Padding = 'tight';

        % Maximize figure window size
        fig.WindowState = 'maximized';    

        % Save figure
        if SAVE_FIGURES
            if FILTER_HIGHPASS_ON
                [~, save_file, ~] = fileparts(filename_stem);
                savefig_multiformat(gcf, SAVE_PATH, [save_file, '_', CHANNEL_KEY{i}, '_blackmanHPfilter', num2str(N_ORDER)])
    %             savefig_multiformat(gcf, SAVE_PATH, [save_file, '_', CHANNEL_KEY{i}, '_movingavgHPfilter', num2str(N_ORDER)])
            else
                [~, save_file, ~] = fileparts(filename_stem);
                savefig_multiformat(gcf, SAVE_PATH, [save_file, '_', CHANNEL_KEY{i}])
            end
        end
    
    end

    
    % Calculate thresholds at each frequency for CAP and ABR
    all_freqs = sort(unique(freq_csv_cache));
    n_freq = numel(all_freqs);
    threshold = zeros(n_freq, 1);
    metric = zeros(NUM_LEVELS, n_freq);
    y_avg_cache_descending = cell(n_freq, 1);
    y_avg_cm_cache_descending = cell(n_freq, 1);
    for j = 1:n_freq
        this_freq = all_freqs(j);
        is_freq = freq_csv_cache==this_freq;
        [A_descending, I] = sort(A_csv_cache(is_freq), 'descend');
        M = cell2mat(y_avg_cache); 
        this_M = M(is_freq, :);% Matrix of size (n_traces x n_samples)
        M_sorted = this_M(I, :);

        [this_threshold, this_metric] = get_thresh_averageABR_oghalai(M_sorted, A_descending, false);
        
        % Cache
        threshold(j) = this_threshold;
        metric(:, j) = this_metric;
        y_avg_cache_descending{j} = M_sorted; % each entry is sorted stack of average traces for a frequency
        
        % Cache microphonic
        M_cm = cell2mat(y_avg_cm_cache); 
        this_M_cm = M_cm(is_freq, :);% Matrix of size (n_traces x n_samples)
        M_cm_sorted = this_M_cm(I, :);
        y_avg_cm_cache_descending{j} = M_cm_sorted;        
    end
    
%     % Reshape cache of single trial data
    X_csv_cache = reshape(X_csv_cache, [NUM_LEVELS, n_freq]);
    wave1peak_cache = reshape(wave1peak_cache, [NUM_LEVELS, n_freq]);
    wave1trough_cache = reshape(wave1trough_cache, [NUM_LEVELS, n_freq]);
    wave1amp_cache = reshape(wave1amp_cache, [NUM_LEVELS, n_freq]);
    wave1lat_peak_cache = reshape(wave1lat_peak_cache, [NUM_LEVELS, n_freq]);
    wave1lat_trough_cache = reshape(wave1lat_trough_cache, [NUM_LEVELS, n_freq]);
    dist_wave1amp = reshape(dist_wave1amp, [NUM_LEVELS, n_freq]);
    wave1std = reshape(wave1std, [NUM_LEVELS, n_freq]);
    
    % Calculate threshold using single trial inner products
    [threshold_D, metric_D] = get_thresh_CAPABR_D_allfreq(X_csv_cache, A_descending, all_freqs, false);
%     [threshold_D, metric_D] = get_thresh_CAPABR_D_allfreq(X_csv_cache, A_descending, all_freqs, true); % Show plots
    
    % Save average traces and their thresholds at each frequency for CAP and ABR
    [~, save_file, ~] = fileparts(filename_stem);
    save_filename = [save_file, '_', CHANNEL_KEY{i}, '_results.mat'];
    save(fullfile(SAVE_PATH, save_filename), 'y_avg_cache_descending', ...
        'y_avg_cm_cache_descending', 'all_freqs', 'A_descending', ...
        'threshold', 'metric', 'threshold_D', 'metric_D', ... % don't save 'X_csv_cache' because too large
        'amp_cm_8khz90db', ...
        'wave1peak_cache', ...
        'wave1trough_cache', ...
        'wave1amp_cache', ...
        'wave1lat_peak_cache', ...
        'wave1lat_trough_cache', ...
        'dist_wave1amp', ...
        'wave1std' ...
        )
    
end


disp('Done')

end