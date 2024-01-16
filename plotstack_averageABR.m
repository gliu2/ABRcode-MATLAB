function outTable = plotstack_averageABR(varargin)
% Plot a stack of average ABRs in CSV file, one stack per frequency. Saves
% plots in multiple image formats.
% Analyzes wave 1 amplitude and latency in highest stimulus waveform.
% Estimate ABR threshold.
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
% 11/17/21 - Added RMS as another output metric for quality control of
% average ABR data.
%
% 10/6/2021 George Liu
% Dependencies: import_averageABRcsv.m, savefig_multiformat.m,
% get_wave1_averageABR.m, get_roc_innerprod_ABR,
% get_thresh_averageABR_liberman.m, get_thresh_averageABR_oghalai.m

%% Constants
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures
SAMPLE_PERIOD = 40.96; % sample period of 40.96 us/sample in average trace
SAMPLE_PERIOD_MS = SAMPLE_PERIOD/1000;
TIME_WINDOW = [1, 2]; % time in ms to analyze portion of ABR trace with inner product method
TIME_BORDER = 0.1; % ms before peak and after trough to include in time window
TIME_BEFOREPEAK = 0.2; % start time window this many ms before wave 1 peak - added 12/6/21
TIME_WINDOW_LENGTH = 1.2; % All time windows are this long (ms) from constant time before wave 1 peak

PLOT_FIGURES = 0; % Hard code if plot figures or not

%% Load average trace data
if nargin == 2
    path = varargin{1};
    filename = varargin{2};
elseif nargin == 0
    [filename,path] = uigetfile('*.csv');
else
    disp('Warning: number of input arguments is not 0 or 2!')
    [filename,path] = uigetfile('*.csv');
end

disp(['Opening ', fullfile(path, filename)])
[M, A_csv, freq_csv] = import_averageABRcsv(filename, path);
num_samples = size(M, 2);
X = SAMPLE_PERIOD_MS * (1:num_samples); % time in ms
A_levels = unique(A_csv);
num_levels = length(A_levels);

% Analyze ABRs for each frequency separately
[freq_unique, ~, ic] = unique(freq_csv);
n_freq = length(freq_unique);

% Initialize cache variables
Frequency = freq_unique;
Stimulus_amplitude = zeros(n_freq, num_levels);
Threshold_liberman = zeros(n_freq, 1);
Threshold_oghalai = zeros(n_freq, 1);
Threshold_innerprod = zeros(n_freq, 1);
Threshold_wave1amp = zeros(n_freq, 1);
Metric_liberman = zeros(n_freq, num_levels);
Metric_oghalai = zeros(n_freq, num_levels);
Metric_innerprod = zeros(n_freq, num_levels);
Metric_innerprod_auc = zeros(n_freq, num_levels);
Metric_wave1amp = zeros(n_freq, num_levels);
Metric_wave1lat = zeros(n_freq, num_levels);
time_window_aroundwave1 = zeros(n_freq*num_levels, 4);
rms_averagetrace = zeros(n_freq, num_levels);

% Obtain y axis scale
ylim_max = [-1200, 1200];
% Check to make sure bounds of plot are within y limits
ymax_data = max(M(:));
ymin_data = min(M(:));
if ymax_data > ylim_max(2)
    ylim_max(2) = ymax_data;
end
if ymin_data < ylim_max(1)
    ylim_max(1) = ymin_data;
end

% Plot ABRs in column stacks for each frequencies, in one tiled layout
if PLOT_FIGURES
    fig = figure;
    t = tiledlayout(num_levels, n_freq, 'TileIndexing', 'columnmajor');
end
for ff = 1:n_freq
    this_freq = freq_unique(ff);
%     disp(['Working on ', num2str(this_freq), ' Hz...'])
    
    % amplitudes and DPOAE values for this frequency only
    A_csv2 = A_csv(ic==ff, :);
    M2 = M(ic==ff, :);
    
%     num_levels = size(M2, 1);
    dB = cell(num_levels, 1);
    dB(:) = {'dB (nV)'};

    [A_descending, I] = sort(A_csv2, 'descend');
    M_sorted = M2(I, :);
    A_descending_cell = cellstr(num2str(A_descending));
    newYlabels=cellfun(@(x,y) [x  ' ' y], A_descending_cell, dB, 'un', 0);

    %% Plot ABRs in vertical stack
%     figure
%     t = tiledlayout(num_levels, 1);
%     s = stackedplot(X, M_sorted', 'Title', plot_title, 'DisplayLabels', newYlabels, 'LineWidth', 1.5);
    
    % Calculate threshold
    [thresh_cache_liberman, this_metric_liberman] = get_thresh_averageABR_liberman(M_sorted, A_descending);
    [thresh_cache_oghalai, this_metric_oghalai] = get_thresh_averageABR_oghalai(M_sorted, A_descending, false);

    % initialize variables
    is_abovenoise = zeros(num_levels, 1);
    for i = 1:num_levels
        y = M_sorted(i, :);
        
        if PLOT_FIGURES
            nexttile(t)
            plot(X, y, 'LineWidth', 3)
            ylim(ylim_max)

            % show ylabels for first column only
            if ff==1
                ylabel(newYlabels{i}, 'FontSize', 24)
                % rotate y label to make horizontal
                hYLabel = get(gca,'YLabel');
                set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right', 'FontSize', 24)
            end
        end

        % Mark wave 1 peak and following trough
        if i==1
            % Get wave 1 peak and following trough at highest stimulus
            % level
            [peak_pt, trough_pt, amp, lat_peak, lat_trough] = get_wave1_averageABR(X, y);
        else
            % Get wave 1 peak and following trough at lower stimulus
            % level. Ensure peak latency does not decrease.
            [peak_pt, trough_pt, amp, lat_peak, lat_trough] = get_wave1_averageABR(X, y, peak_pt(1), trough_pt(1));
        end
        
        if PLOT_FIGURES
            hold on
            CIRCLE_SIZE = 54; % default 36
            scatter(peak_pt(1), peak_pt(2), CIRCLE_SIZE, 'green', 'filled') % peak
            scatter(trough_pt(1), trough_pt(2), CIRCLE_SIZE, 'red', 'filled') % trough
            hold off

            % Display wave 1 measurements in corner of plot
            xL=xlim;
            yL=ylim;
            str = sprintf('P-P_I = %.0f nV', amp);
            if A_descending(i) > thresh_cache_liberman
                text(0.99*xL(2), 0.99*yL(1), str,'HorizontalAlignment','right','VerticalAlignment','bottom', 'FontWeight', 'bold')
                is_abovenoise(i,1)=1;
            else
                text(0.99*xL(2), 0.99*yL(1), str,'HorizontalAlignment','right','VerticalAlignment','bottom')
                is_abovenoise(i,1)=0;
            end

            % Remove extraneous axis tick marks and x-axis from all but bottom
            % tile
            set(gca,'box','off')
            if i~=num_levels
                set(gca,'xtick',[])
                set(gca,'xticklabel',[])
                set(gca,'XColor','none')
            end

            % Remove y-axis labels from all but first column
            if ff~=1
                set(gca,'yticklabel',[])
            end

            % Set text size of current axis
            set(gca,'FontSize',24)

            % Show frequency as title above top-most tile of each column
            if i==1
                title([num2str(this_freq), ' Hz'])
            end
        end
        
        % Cache wave 1 amplitude and latency for highest stimulus level
        Metric_wave1amp(ff, i) = amp; % nV
        Metric_wave1lat(ff, i) = lat_peak; % ms
        % Cache wave 1 peak and trough latency to estimate time window for
        % calculating inner product around wave 1
        count = (ff - 1)*num_levels + i;
        time_window_aroundwave1(count, 1) = lat_peak; % ms
        time_window_aroundwave1(count, 2) = lat_trough; % ms
        time_window_aroundwave1(count, 3) = this_freq;
        time_window_aroundwave1(count, 4) = A_descending(i);
        
    end

    xlabel('Time (ms)', 'FontSize', 24)
    
    % Calculate threshold based on wave 1 amplitude
    noise_std = std(M_sorted(end, :));
    wave1_thresh = get_thresh_wave1amp_oghalai(Metric_wave1amp(ff, :), A_descending, noise_std, false);
    
    % Cache variables
    Threshold_wave1amp(ff) = wave1_thresh;
    Stimulus_amplitude(ff, :) = A_descending';
    Threshold_liberman(ff) = thresh_cache_liberman; % dB
    Threshold_oghalai(ff) = thresh_cache_oghalai;
    Metric_liberman(ff, :) = this_metric_liberman;
    Metric_oghalai(ff, :) = this_metric_oghalai;
    
    % Calculate RMS of average trace as quality control metric
    rms_averagetrace(ff, :) = rms(M_sorted, 2);
    
end


if PLOT_FIGURES
    t.TileSpacing = 'none';
    t.Padding = 'tight';
    %     plot_title = ['ABR @ ', num2str(this_freq), ' Hz'];
    %     title(t, plot_title)
    %     ylabel(t, 'Amplitude (nV)')


    % Maximize figure window size
    fig.WindowState = 'maximized';

    % Save figure
    %     disp('Saving figure')
    [~, save_file, ~] = fileparts(filename);
    savefig_multiformat(gcf, SAVE_PATH, save_file)
end

% Obtain wave 1 time window (1-2 ms) single trace inner product threshold and AUC level functions
% Edit time window to ensure that it is in units of ms, and that a border
% is placed around to capture some data before peak and after trough
% time_window_aroundwave1(:,1) = time_window_aroundwave1(:,1)*SAMPLE_PERIOD_MS - TIME_BORDER; % window start time (ms)
% time_window_aroundwave1(:,2) = time_window_aroundwave1(:,2)*SAMPLE_PERIOD_MS + TIME_BORDER; % window end time (ms)
% time_window_aroundwave1(:,1) = time_window_aroundwave1(:,1)*SAMPLE_PERIOD_MS - TIME_BEFOREPEAK; % window start time (ms)
% if time_window_aroundwave1(:,1) < 0
%     time_window_aroundwave1(:,1) = 0;
% end
% time_window_aroundwave1(:,2) = time_window_aroundwave1(:,1) + TIME_WINDOW_LENGTH; % window end time (ms)
% time_window_aroundwave1(:,1:2) = time_window_aroundwave1(:,1:2)*SAMPLE_PERIOD_MS; % wave 1 peak and trough latencies are already in units of time (ms)
% [auc_window_cache, ip_window_mean_cache, ip_window_ste_cache, thresh_window_cache_innerprod, thresh_window_criterion, ~, ~, ~, ~] = get_roc_innerprod_ABR(fullfile(path, filename), false, time_window_aroundwave1);
signal_metrics_thresholds = get_roc_innerprod_ABR(fullfile(path, filename), false, time_window_aroundwave1); % Don't plot debugging figures
% signal_metrics_thresholds = get_roc_innerprod_ABR(fullfile(path, filename), true, time_window_aroundwave1); % Plot debugging figures
if length(signal_metrics_thresholds)==2
    [auc_window_cache, ip_window_mean_cache, ip_window_ste_cache, thresh_window_cache_innerprod, thresh_window_criterion, ~, ~, ~, ~, dist_wave1amp_window_cache, dist_innerprod_window_cache, auc_cache_wave1] = signal_metrics_thresholds{2}{:};
    Threshold_innerprod_window = thresh_window_cache_innerprod;
    Threshold_innerprod_window_criterion = thresh_window_criterion;
    Metric_innerprod_window = ip_window_mean_cache;
    Metric_innerprod_window_ste = ip_window_ste_cache;
    Metric_innerprod_window_auc = auc_window_cache;
    Distribution_innerprod_window = dist_innerprod_window_cache;
end

% Obtain single trace inner product threshold and AUC level functions
% [auc_cache, ip_mean_cache, ip_ste_cache, thresh_cache_innerprod, thresh_criterion, d_prime_cache_rms, thresh_d_prime_rms, d_prime_cache_z, thresh_d_prime_z] = get_roc_innerprod_ABR(fullfile(path, filename));
[auc_cache, ip_mean_cache, ip_ste_cache, thresh_cache_innerprod, thresh_criterion, d_prime_cache_rms, thresh_d_prime_rms, d_prime_cache_z, thresh_d_prime_z, dist_wave1amp_cache, dist_innerprod_cache, auc_cache_wave1] = signal_metrics_thresholds{1}{:};
Threshold_innerprod = thresh_cache_innerprod;
Threshold_innerprod_criterion = thresh_criterion; % 11/12/21 - Assign inner product threshold calculated with criterion to inner product label, and Wilcoxin rank-sum IP threshold to AUC inner product label
Metric_innerprod = ip_mean_cache;
Metric_innerprod_ste = ip_ste_cache;
Metric_innerprod_auc = auc_cache;
Metric_auc_wave1 = auc_cache_wave1;
Distribution_innerprod = dist_innerprod_cache;
Distribution_wave1amp = dist_wave1amp_cache;
Threshold_d_rms = thresh_d_prime_rms; % 12/6/21 - added d' signal detection sensitivity metrics
Threshold_d_z = thresh_d_prime_z;
Metric_d_rms = d_prime_cache_rms;
Metric_d_z = d_prime_cache_z;


% Write cached variables to excel table
% Structure of table is: filename, frequency, method, threshold, metric
% level function.
% Listing of filenames
[~, baseFileNameNoExt, ~] = fileparts(filename);
Filenames = cell(n_freq, 1);
Filenames(:) = {baseFileNameNoExt};  
Filenames = convertCharsToStrings(Filenames);

% Metric names
Metric_names_amplitude = repmat("Amplitude", n_freq, 1);
Metric_names_liberman = repmat("Liberman", n_freq, 1);
Metric_names_oghalai = repmat("Oghalai", n_freq, 1);
Metric_names_innerprod = repmat("Innerprod", n_freq, 1);
Metric_names_innerprod_auc = repmat("Innerprod_auc", n_freq, 1);
Metric_names_innerprod_window = repmat("Innerprod_window", n_freq, 1);
Metric_names_innerprod_window_auc = repmat("Innerprod_window_auc", n_freq, 1);
Metric_names_wave1amp = repmat("Wave1amp", n_freq, 1);
Metric_names_wave1amp_auc = repmat("Wave1amp_auc", n_freq, 1);
Metric_names_wave1lat = repmat("Wave1lat", n_freq, 1);
Metric_names_rms = repmat("rms", n_freq, 1);
Metric_names_d_rms = repmat("D_rms", n_freq, 1);
Metric_names_d_z = repmat("D_z", n_freq, 1);
blank_space = nan(n_freq, 1);
blank_space2 = cell(n_freq, 1); 
blank_space2(:) = {NaN};

% Stimulus amplitudes in place of metrics
outTable_amplitudes = table(Filenames, Frequency, Metric_names_amplitude, blank_space , blank_space2,...
    Stimulus_amplitude);
outTable_amplitudes.Properties.VariableNames{'Metric_names_amplitude'} = 'Metric';
outTable_amplitudes.Properties.VariableNames{'blank_space'} = 'Threshold';
outTable_amplitudes.Properties.VariableNames{'blank_space2'} = 'Distribution';
outTable_amplitudes.Properties.VariableNames{'Stimulus_amplitude'} = 'Metric_values';
% Liberman method results
outTable_liberman = table(Filenames, Frequency, Metric_names_liberman, Threshold_liberman, blank_space2 ,...
    Metric_liberman);
outTable_liberman.Properties.VariableNames{'Metric_names_liberman'} = 'Metric';
outTable_liberman.Properties.VariableNames{'Threshold_liberman'} = 'Threshold';
outTable_liberman.Properties.VariableNames{'blank_space2'} = 'Distribution';
outTable_liberman.Properties.VariableNames{'Metric_liberman'} = 'Metric_values';
% Oghalai method results
outTable_oghalai = table(Filenames, Frequency, Metric_names_oghalai, Threshold_oghalai, blank_space2 ,...
    Metric_oghalai);
outTable_oghalai.Properties.VariableNames{'Metric_names_oghalai'} = 'Metric';
outTable_oghalai.Properties.VariableNames{'Threshold_oghalai'} = 'Threshold';
outTable_oghalai.Properties.VariableNames{'blank_space2'} = 'Distribution';
outTable_oghalai.Properties.VariableNames{'Metric_oghalai'} = 'Metric_values';
% Wave 1 latency
outTable_wave1lat = table(Filenames, Frequency, Metric_names_wave1lat, Threshold_wave1amp, blank_space2 ,...
    Metric_wave1lat);
outTable_wave1lat.Properties.VariableNames{'Metric_names_wave1lat'} = 'Metric';
outTable_wave1lat.Properties.VariableNames{'Threshold_wave1amp'} = 'Threshold';
outTable_wave1lat.Properties.VariableNames{'blank_space2'} = 'Distribution';
outTable_wave1lat.Properties.VariableNames{'Metric_wave1lat'} = 'Metric_values';
% RMS of 30 dB average trace
outTable_rms = table(Filenames, Frequency, Metric_names_rms, Threshold_wave1amp, blank_space2 ,...
    rms_averagetrace);
outTable_rms.Properties.VariableNames{'Metric_names_rms'} = 'Metric';
outTable_rms.Properties.VariableNames{'Threshold_wave1amp'} = 'Threshold';
outTable_rms.Properties.VariableNames{'blank_space2'} = 'Distribution';
outTable_rms.Properties.VariableNames{'rms_averagetrace'} = 'Metric_values';

if ~isempty(Threshold_innerprod)
    % Wave 1 amplitude method results
    outTable_wave1amp = table(Filenames, Frequency, Metric_names_wave1amp, Threshold_wave1amp,...
        Distribution_wave1amp, Metric_wave1amp);
    outTable_wave1amp.Properties.VariableNames{'Metric_names_wave1amp'} = 'Metric';
    outTable_wave1amp.Properties.VariableNames{'Threshold_wave1amp'} = 'Threshold';
    outTable_wave1amp.Properties.VariableNames{'Distribution_wave1amp'} = 'Distribution';
    outTable_wave1amp.Properties.VariableNames{'Metric_wave1amp'} = 'Metric_values';
    % Wave 1 amplitude AUC results
    outTable_wave1amp_AUC = table(Filenames, Frequency, Metric_names_wave1amp_auc, Threshold_wave1amp,...
        Distribution_wave1amp, Metric_auc_wave1);
    outTable_wave1amp_AUC.Properties.VariableNames{'Metric_names_wave1amp_auc'} = 'Metric';
    outTable_wave1amp_AUC.Properties.VariableNames{'Threshold_wave1amp'} = 'Threshold';
    outTable_wave1amp_AUC.Properties.VariableNames{'Distribution_wave1amp'} = 'Distribution';
    outTable_wave1amp_AUC.Properties.VariableNames{'Metric_auc_wave1'} = 'Metric_values';
    % Inner product method results 
%     outTable_innerprod = table(Filenames, Frequency, Metric_names_innerprod, Threshold_innerprod, ...
%         Metric_innerprod);
    outTable_innerprod = table(Filenames, Frequency, Metric_names_innerprod, Threshold_innerprod_criterion, ... % 11/12/21 - Assign inner product threshold calculated with criterion to inner product label, and Wilcoxin rank-sum IP threshold to AUC inner product label
        Distribution_innerprod, Metric_innerprod);
    outTable_innerprod.Properties.VariableNames{'Metric_names_innerprod'} = 'Metric';
    outTable_innerprod.Properties.VariableNames{'Threshold_innerprod_criterion'} = 'Threshold';
    outTable_innerprod.Properties.VariableNames{'Distribution_innerprod'} = 'Distribution';
    outTable_innerprod.Properties.VariableNames{'Metric_innerprod'} = 'Metric_values';
    % Inner product AUC method results
    outTable_innerprod_AUC = table(Filenames, Frequency, Metric_names_innerprod_auc, Threshold_innerprod, ...
        Distribution_innerprod, Metric_innerprod_auc);
    outTable_innerprod_AUC.Properties.VariableNames{'Metric_names_innerprod_auc'} = 'Metric';
    outTable_innerprod_AUC.Properties.VariableNames{'Threshold_innerprod'} = 'Threshold';
    outTable_innerprod_AUC.Properties.VariableNames{'Distribution_innerprod'} = 'Distribution';
    outTable_innerprod_AUC.Properties.VariableNames{'Metric_innerprod_auc'} = 'Metric_values';
    % Inner product time window method results
    outTable_innerprod_window = table(Filenames, Frequency, Metric_names_innerprod_window, Threshold_innerprod_window_criterion, ...
        Distribution_innerprod_window, Metric_innerprod_window);
    outTable_innerprod_window.Properties.VariableNames{'Metric_names_innerprod_window'} = 'Metric';
    outTable_innerprod_window.Properties.VariableNames{'Threshold_innerprod_window_criterion'} = 'Threshold';
    outTable_innerprod_window.Properties.VariableNames{'Distribution_innerprod_window'} = 'Distribution';
    outTable_innerprod_window.Properties.VariableNames{'Metric_innerprod_window'} = 'Metric_values';
    % Inner product time window AUC results
    outTable_innerprod_window_AUC = table(Filenames, Frequency, Metric_names_innerprod_window_auc, Threshold_innerprod_window, ...
        Distribution_innerprod_window, Metric_innerprod_window_auc);
    outTable_innerprod_window_AUC.Properties.VariableNames{'Metric_names_innerprod_window_auc'} = 'Metric';
    outTable_innerprod_window_AUC.Properties.VariableNames{'Threshold_innerprod_window'} = 'Threshold';
    outTable_innerprod_window_AUC.Properties.VariableNames{'Distribution_innerprod_window'} = 'Distribution';
    outTable_innerprod_window_AUC.Properties.VariableNames{'Metric_innerprod_window_auc'} = 'Metric_values';
    % D prime 'RMS' : d' = u_signal - u_noise / rms(std_signal, std_noise) - 
    outTable_d_rms = table(Filenames, Frequency, Metric_names_d_rms, Threshold_d_rms, ...
        blank_space2, Metric_d_rms);
    outTable_d_rms.Properties.VariableNames{'Metric_names_d_rms'} = 'Metric';
    outTable_d_rms.Properties.VariableNames{'Threshold_d_rms'} = 'Threshold';
    outTable_d_rms.Properties.VariableNames{'blank_space2'} = 'Distribution';
    outTable_d_rms.Properties.VariableNames{'Metric_d_rms'} = 'Metric_values';
    % D prime 'z' : d' = z_noise(criterion) - z_signal(criterion)  
    outTable_d_z = table(Filenames, Frequency, Metric_names_d_z, Threshold_d_z, ...
        blank_space2, Metric_d_z);
    outTable_d_z.Properties.VariableNames{'Metric_names_d_z'} = 'Metric';
    outTable_d_z.Properties.VariableNames{'Threshold_d_z'} = 'Threshold';
    outTable_d_z.Properties.VariableNames{'blank_space2'} = 'Distribution';
    outTable_d_z.Properties.VariableNames{'Metric_d_z'} = 'Metric_values';
else
    outTable_innerprod = [];
    outTable_innerprod_AUC = [];
    outTable_innerprod_window = [];
    outTable_innerprod_window_AUC = [];
    outTable_wave1amp_AUC = [];
    outTable_noiseamp = [];
    outTable_d_rms = [];
    outTable_d_z = [];
    % Wave 1 amplitude method results
    outTable_wave1amp = table(Filenames, Frequency, Metric_names_wave1amp, Threshold_wave1amp, blank_space2 ,...
    Metric_wave1amp);
    outTable_wave1amp.Properties.VariableNames{'Metric_names_wave1amp'} = 'Metric';
    outTable_wave1amp.Properties.VariableNames{'Threshold_wave1amp'} = 'Threshold';
    outTable_wave1amp.Properties.VariableNames{'blank_space2'} = 'Distribution';
    outTable_wave1amp.Properties.VariableNames{'Metric_wave1amp'} = 'Metric_values';
end

outTable = [outTable_amplitudes; outTable_liberman; outTable_oghalai; outTable_innerprod; outTable_innerprod_AUC; ...
    outTable_innerprod_window; outTable_innerprod_window_AUC; outTable_wave1amp; outTable_wave1amp_AUC; outTable_wave1lat; outTable_rms; outTable_d_rms; outTable_d_z];

% % Save excel sheet output
% table_filename = [baseFileNameNoExt, '_output.xlsx'];
% writetable(outTable, fullfile(SAVE_PATH, table_filename), 'Sheet', 1, 'Range', 'A1')

% disp('Done')

end