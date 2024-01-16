% function plot_averageABR_prepost_RMSsingletraces(varargin)
% Plot average ABR waveform (+/- pre and post noise exposure) at a single
% stimulus and frequency in CSV file, with second plot of RMS of single traces to evaluate for artifacts. 
%
% Default is for user to select file in dialog box when no input parameters
% are specified.
%
% Optional input syntax for specifying input file: plot_averageABR_RMSsingletraces(path, filename)
%
% 9/27/2022 George Liu
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

XLIM_RMS = [0, 9000];
YLIM_RMS = [0, 500];

HISTOGRAM_BINS = 50;
N_SINGLETRACES = 40;
SCATTER_CIRCLE_SIZE = 54; % default 36
FONTSIZE = 14;

PLOT_REJECTION_SINGLETRACES = 0;
SAVE_FIGURES = 0;

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

% Initialize cache for single trace RMS values for later plotting of RMS vs
% single trace index
singletraces_rms_cache = cell(num_mice, num_dates_per_mouse);
singletraces_xcorr_cache = cell(num_mice, num_dates_per_mouse);

% initialize cache for storing averaged ABR waveforms after artifact
% rejection based on threshold of RMS criterion 
averageABR = cell(num_mice, num_dates_per_mouse);
averageABR_artifactrejection = cell(num_mice, num_dates_per_mouse);
averageABR_artifact_filtered = cell(num_mice, num_dates_per_mouse);
averageABR_noise = cell(num_mice, num_dates_per_mouse);

% Initialize variables for plotting correlation of noise ABR RMS with wave
% 1 amp in stimulus case
rms_no_stimulus_averageABR = zeros(num_mice, num_dates_per_mouse);
wave1amp = zeros(num_mice, num_dates_per_mouse);

% % Initialize and visualize Butterworth filter for processing average ABR
% % trace, using Butterworth filter parameters from Liberman group.
% LOWFREQ = 200; % Hz
% HIGHFREQ = 10000; % Hz
% N_ORDER = 2;
% d = designfilt('bandpassiir','FilterOrder', N_ORDER, ...
% 'HalfPowerFrequency1', LOWFREQ,'HalfPowerFrequency2', HIGHFREQ, ...
% 'SampleRate', SAMPLING_RATE_HZ);
% fvt = fvtool(d, 'Fs', SAMPLING_RATE_HZ);
% legend(fvt,'bandpassiir')

% 10-4-22: as an alternative, use Blackman filter
% NUM_MS_WAVE1_PEAK = 0.47;
NUM_MS_WAVE1_PEAK = 0.3;
L = round(SAMPLING_RATE*NUM_MS_WAVE1_PEAK);
w = blackman(L);
w = w/sum(w); % Normalization

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

        averageABR{j, i} = M2; % store average ABR trace
        
        % 10-11-2022: Obtain average ABR in noise at lowest stimlus level
        lowest_stimulus_level = min(A_csv);
        is_lowest_stimulus = A_csv==lowest_stimulus_level;
        averageABR_noise{j, i} = M(is_freq & is_lowest_stimulus, :);
        
        % 9-21-2022: Obtain single ABR traces at the specified frequency
        % and intensity
        input_pathfile = fullfile(this_path, filename);
        [this_Xcsv, x_ms] = import_singleABR_frequency(input_pathfile, THIS_FREQ, THIS_INTENSITY); % SAMPLES x m/2 matrix
        
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
        
        % Calculate RMS of each single trace
        singletraces_rms = rms(this_Xcsv, 1); % calculate RMS for each column (single trace) 
        
        % Cache RMS of single traces
        singletraces_rms_cache{j, i} = singletraces_rms;
        
        % Calculate cross-correlation of single trace with average trace,
        % to determine large driven responses in single trace "artifacts"
        % (that are not actually artifacts since shape of response is
        % similar to average ABR)
        MAX_LAG = 20;
        num_traces = numel(singletraces_rms);
        this_xcorr = zeros(num_traces , 1);
        for k=1:num_traces
            this_trace_xcorr = xcorr(this_Xcsv(:, k), mean(this_Xcsv, 2), MAX_LAG, 'coeff');
            this_xcorr(k) = max(this_trace_xcorr(:));
        end
        singletraces_xcorr_cache{j, i} = this_xcorr;
        
        % Determine RMS cutoff for artifact rejection using mean + std
        % criteria
        median_singletraces_rms = median(singletraces_rms); 
        std_singletraces_rms = std(singletraces_rms);
        cutoff = median_singletraces_rms + N_STD_CUTOFF*std_singletraces_rms;

        % 9-28-22: Determine RMS cutoff for artifact rejection using measure of
        % RMS in average of rejected single traces
        if PLOT_REJECTION_SINGLETRACES
            artifact_rms = 0;
            count = 0;
    %         ARTIFACT_RMS_MINIMUM = mean_singletraces_rms - std_singletraces_rms;
    %         while artifact_rms < ARTIFACT_RMS_MINIMUM

    %         this_ylim = [min(this_Xcsv(:)), max(this_Xcsv(:))];
            this_ylim = 1.5*[-10^4, 10^4];
            fig2 = figure;
            t2 = tiledlayout(NUM_ITERATIONS, 3, 'TileIndexing', 'rowmajor');
            while count < NUM_ITERATIONS
                count = count + 1;
                disp(['Iterating on RMS cutoff, attempt ', num2str(count)])

                % Update cutoff for rejecting single traces
                num_std = N_STD_CUTOFF - count + 1;
                cutoff = median_singletraces_rms + num_std*std_singletraces_rms;

                isbelow_cutoff = singletraces_rms < cutoff;


                % plot RMS cutoff on histogram plot, average of non-rejected
                % traces, and average of rejected traces


                % Histogram
                nexttile(t2)
                histogram(singletraces_rms, HISTOGRAM_BINS)
                xlabel('RMS (nV)')
                ylabel('Frequency')
                hold on
                xline(median_singletraces_rms, '--', 'median')
                xline(cutoff, '--', ['cutoff = median + ', num2str(num_std), '*std'])
                if count==1
                    title('RMS histogram')
                end
                % 10-11-22: Draw red vertical line at RMS of noise
                % average ABR trace
                rms_no_stimulus_averageABR(j, i) = rms(averageABR_noise{j, i});
                xline(rms_no_stimulus_averageABR(j, i), '--r', 'RMS of average ABR with no stimulus')
                
    %             % Make uniform x and y limits
    %             xlim(XLIM_RMS)
    %             ylim(YLIM_RMS)

                % Set text size of current axis
                set(gca,'FontSize', FONTSIZE)


                % Average of non-rejected traces
                averageABR_removed_artifacts = mean(this_Xcsv(:, isbelow_cutoff), 2);
                ind_isbelow_cutoff = find(isbelow_cutoff);
                nexttile(t2)

                % Plot single traces
                n_traces = sum(isbelow_cutoff);
                cmap = summer(n_traces);
                hold on
                for k=1:n_traces
                    plot(x_ms, this_Xcsv(:, ind_isbelow_cutoff(k)), 'Color', [cmap(k, 1), cmap(k, 2), cmap(k, 3), 0.8]) % transparency set to 4th color value
                end

                % Plot raw and non-rejected traces
                p1=plot(X, averageABR{j, i}, 'blue', 'LineWidth', 3);
                p2=plot(x_ms, averageABR_removed_artifacts, 'magenta', 'LineWidth', 3);

                if count==1
                    legend([p1, p2], {'Raw average', 'Artifact rejection'}, 'location', 'southeast', 'FontSize', FONTSIZE)
                    title('Artifact rejection')
                    colorbar
                end

                xlabel('Time (ms)')
                ylabel('Amplitude (nV)')
                ylim(this_ylim)
                set(gca,'FontSize', FONTSIZE)



                % Average of rejected traces
                artifact_average = mean(this_Xcsv(:, ~isbelow_cutoff), 2);
                ind_isabove_cutoff = find(~isbelow_cutoff);
                artifact_rms = rms(artifact_average);
                nexttile(t2)
                % Plot single traces
                n_traces = sum(~isbelow_cutoff);
                cmap = summer(n_traces);
                hold on
                for k=1:n_traces
                    plot(x_ms, this_Xcsv(:, ind_isabove_cutoff(k)), 'Color', [cmap(k, 1), cmap(k, 2), cmap(k, 3), 0.8]) % transparency set to 4th color value
                end

                p1=plot(X, averageABR{j, i}, 'blue', 'LineWidth', 3);
                hold on
                p2=plot(x_ms, artifact_average, 'magenta', 'LineWidth', 3);

                xlabel('Time (ms)')
                ylabel('Amplitude (nV)')
                if count==1
                    legend([p1, p2], {'Raw average', 'Rejected average'}, 'location', 'southeast', 'FontSize',16)
                    title(['Rejected average'])
                end
                ylim(this_ylim)
                set(gca,'FontSize', FONTSIZE)

            end

            t2.TileSpacing = 'none';
            % t.TileSpacing = 'tight';
            t2.Padding = 'tight';

            % Maximize figure window size
            fig2.WindowState = 'maximized';
            
            % Save figure of 3-column tiled plot of (1) histogram with RMS
            % cutoff, (2) average trace after artifact rejection, and (3)
            % average of rejected traces
            if SAVE_FIGURES
                save_file = [this_mouse_label, '_histogram_rejectartifacts_rejectedaverage_', this_date];
                savefig_multiformat(gcf, SAVE_PATH, save_file)
            end
        end
        
        % Store average ABR trace with artifact single traces rejected
        isbelow_cutoff = singletraces_rms < cutoff;
        averageABR_artifactrejection{j, i} = mean(this_Xcsv(:, isbelow_cutoff), 2);
        
%         % Store average ABR trace after artifact rejection, AND bandpass
%         % butterworth filtering of average trace, using zero-phase 
%         % first-order Butterworth filter with a pass band from 0.2 to 10 kHz
%         averageABR_artifact_filtered{j, i} = filter(d, averageABR_artifactrejection{j, i});

        % Apply Blackman filter to average trace
        averageABR_artifact_filtered{j, i} = conv(averageABR_artifactrejection{j, i}, w, 'same');
        
        
        % Plot histogram of single trace RMS value distribution
        nexttile(t)
        histogram(singletraces_rms, HISTOGRAM_BINS)
        hold on
        xline(median_singletraces_rms, '--', 'median')
        xline(cutoff, '--', ['median + ', num2str(N_STD_CUTOFF), '*std'])

        % 10-11-22: Draw red vertical line at RMS of noise
        % average ABR trace
        rms_no_stimulus_averageABR(j, i) = rms(averageABR_noise{j, i});
        xline(rms_no_stimulus_averageABR(j, i), '--r', 'RMS of average ABR with no stimulus')
        
        % Make uniform x and y limits
        xlim(XLIM_RMS)
        ylim(YLIM_RMS)
        
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
ylabel(t, 'Frequency (count)', 'FontSize', 24)
xlabel(t, 'RMS (nV)', 'FontSize', 24)

% Maximize figure window size
fig.WindowState = 'maximized';

% Save figure
%     disp('Saving figure')
% save_file = [this_mouse_label, '_overlaysingletraces'];
base_file = '';
for i=1:num_mice
    if i==1
        base_file = these_mice_labels{i};
    else
        base_file = [base_file, '_', these_mice_labels{i}];
    end
end
if SAVE_FIGURES
    save_file = [base_file, '_singletraceRMShistogram'];
    savefig_multiformat(gcf, SAVE_PATH, save_file)
end

%% Plot average ABR waveforms after artifact rejection
% Plot ABRs in one tiled layout
fig = figure;
t = tiledlayout(num_mice, num_dates_per_mouse, 'TileIndexing', 'rowmajor');
for j=1:num_mice
    for i=1:num_dates_per_mouse
        nexttile(t)
        y = averageABR{j, i}; 
        p1 = plot(X, y, 'blue', 'LineWidth', 3);
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
        
        % Overlay plots of average ABR after (1) rejecting single trace
        % artifacts using RMS threshold criterion and (2) subsequent
        % average ABR bandpass filtering
        hold on
        p2 = plot(x_ms, averageABR_artifactrejection{j, i}, 'cyan', 'LineWidth', 3);
        p3 = plot(x_ms, averageABR_artifact_filtered{j, i}, 'magenta', 'LineWidth', 3);
        
        % Show wave 1 peak and trough for original unprocessed average ABR
        % trace
        hold on
        scatter(peak_pt(1), peak_pt(2), SCATTER_CIRCLE_SIZE, 'green', 'filled') % peak
        scatter(trough_pt(1), trough_pt(2), SCATTER_CIRCLE_SIZE, 'red', 'filled') % trough
        hold off
        
        if i==1 && j==1
            legend([p1, p2, p3], {'Raw average', 'Artifact rejection', ['Rejection + Blackman (', num2str(NUM_MS_WAVE1_PEAK), ' ms)']}, 'location', 'southeast', 'FontSize',16)
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
plot_title = ['ABR with artifacts above ', num2str(N_STD_CUTOFF), ' std RMS rejected @ ', num2str(THIS_FREQ), ' Hz, ' num2str(THIS_INTENSITY), ' dB'];
title(t, plot_title, 'FontSize', 24)
ylabel(t, 'Amplitude (nV)', 'FontSize', 24)
xlabel(t, 'Time (ms)', 'FontSize', 24)

% Maximize figure window size
fig.WindowState = 'maximized';

% Save figure
%     disp('Saving figure')
if SAVE_FIGURES
    save_file2 = [base_file, '_averageABR_rejectartifactsRMS', num2str(N_STD_CUTOFF), 'std_blackman', strrep(num2str(NUM_MS_WAVE1_PEAK), '.', 'p')];
    savefig_multiformat(gcf, SAVE_PATH, save_file2)
end

%% 10-4-22: Plot RMS vs index of single trace to assess for systematic variations in artifacts with time under anesthesia 
fig = figure;
t = tiledlayout(num_mice, num_dates_per_mouse, 'TileIndexing', 'rowmajor');

fig2 = figure;
t2 = tiledlayout(num_mice, num_dates_per_mouse, 'TileIndexing', 'rowmajor');

% Determine max RMS for y-limits
max_rms = cellfun(@max, singletraces_rms_cache);
% max_rms_all = max(max_rms(:));
max_rms_all = 8000;
for j=1:num_mice
    for i=1:num_dates_per_mouse
        % Plot RMS vs single trace index
        nexttile(t)
        y_rms = singletraces_rms_cache{j, i}; 
        num_traces = numel(y_rms);
        x_rms = 1:num_traces;
        p1 = scatter(x_rms, y_rms, SCATTER_CIRCLE_SIZE, spring(num_traces), "filled");
        grid on
%         ylim(ylim_max)
        xlim([0, num_traces])
        ylim([0, max_rms_all])

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
                    
        else
            % Remove y-axis labels from all but first column
            set(gca,'yticklabel',[])
        end
        
        
        % Plot RMS vs cross correlation (normalized) with average ABR trace
        % to identify large driven response single traces
        nexttile(t2)
        x_xcorr = singletraces_xcorr_cache{j, i};
        p2 = scatter(x_xcorr, y_rms, SCATTER_CIRCLE_SIZE, spring(num_traces), "filled");
        grid on
        xlim([-1, 1])
        ylim([0, max_rms_all])
        set(gca,'FontSize',24) % Set text size of current axis
        
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
        else
            % Remove y-axis labels from all but first column
            set(gca,'yticklabel',[])
        end

    end
end


% Adjust tiled layout for RMS vs single trace index plot
t.TileSpacing = 'none';
% t.TileSpacing = 'tight';
t.Padding = 'tight';
plot_title = ['RMS per single trace @ ', num2str(THIS_FREQ), ' Hz, ' num2str(THIS_INTENSITY), ' dB, ', num2str(num_traces), ' traces'];
title(t, plot_title, 'FontSize', 24)
ylabel(t, 'RMS (nV)', 'FontSize', 24)
xlabel(t, 'Single trace index', 'FontSize', 24)

% Maximize figure window size
fig.WindowState = 'maximized';

% Save figure
%     disp('Saving figure')
if SAVE_FIGURES
    save_file3 = [base_file, '_RMSpersingletrace'];
    savefig_multiformat(gcf, SAVE_PATH, save_file3)
end



% Adjust tiled layout for RMS vs cross correlation scatter plot
t2.TileSpacing = 'none';
% t.TileSpacing = 'tight';
t2.Padding = 'tight';
plot_title2 = ['RMS vs cross correlation @ ', num2str(THIS_FREQ), ' Hz, ' num2str(THIS_INTENSITY), ' dB, ', num2str(num_traces), ' traces'];
title(t2, plot_title2, 'FontSize', 24)
ylabel(t2, 'RMS (nV)', 'FontSize', 24)
xlabel(t2, 'Cross correlation (normalized)', 'FontSize', 24)

% Maximize figure window size
fig2.WindowState = 'maximized';

% Save figure
%     disp('Saving figure')
if SAVE_FIGURES
    save_file4 = [base_file, '_RMSvsXcorr'];
    savefig_multiformat(gcf, SAVE_PATH, save_file4)
end

%% 10-11-2022 Assess for correlation between noise average ABR RMS and wave 1 amplitude
num_data_points = numel(rms_no_stimulus_averageABR);
% Collect average ABR RMS with stimulus
rms_average_abr = zeros(num_mice, num_dates_per_mouse);
max_amp_noise_abr = zeros(num_mice, num_dates_per_mouse);
group_name = cell(num_data_points, 1);
for j=1:num_mice
    for i=1:num_dates_per_mouse
        ind = i + (j-1)*num_dates_per_mouse;
        rms_average_abr(j, i) = rms(averageABR{j, i});
        group_name{ind, 1} = these_mice_labels{j};
        max_amp_noise_abr(j, i) = max(averageABR_noise{j, i}) - min(averageABR_noise{j, i});
    end
end

grouped_rms_no_stimulus_averageABR = reshape(rms_no_stimulus_averageABR', [num_data_points, 1]);
grouped_rms_average_abr = reshape(rms_average_abr', [num_data_points, 1]);
grouped_wave1amp = reshape(wave1amp', [num_data_points, 1]);
grouped_max_amp_noise_abr = reshape(max_amp_noise_abr', [num_data_points, 1]);

% Scatter plot of wave 1 amp vs noise RMS
fig3 = figure;
t4 = tiledlayout(2,2);
nexttile(t4)
mice_colors = hsv(num_mice);
gscatter(grouped_rms_no_stimulus_averageABR, grouped_wave1amp, group_name, mice_colors, '.', SCATTER_CIRCLE_SIZE/2);
xlabel('RMS of average ABR without stimulus (nV)')
ylabel(['Wave 1 amplitude @ ', num2str(THIS_INTENSITY), ' dB (nV)'])
title(['Wave 1 amplitude vs noise RMS @ ', num2str(THIS_FREQ), ' Hz'])
set(gca,'FontSize', FONTSIZE) % Set text size of current axis

% Scatter plot of average ABR RMS with and without stimulus
nexttile(t4)
gscatter(grouped_rms_no_stimulus_averageABR, grouped_rms_average_abr, group_name, mice_colors, '.', SCATTER_CIRCLE_SIZE/2);
xlabel('RMS of average ABR without stimulus (nV)')
ylabel(['RMS of average ABR @ ', num2str(THIS_INTENSITY), ' dB (nV)'])
title(['RMS with and without stimulus @ ', num2str(THIS_FREQ), ' Hz'])
set(gca,'FontSize', FONTSIZE) % Set text size of current axis

% Scatter plot of wave 1 amp vs noise max amp
nexttile(t4)
gscatter(grouped_max_amp_noise_abr, grouped_wave1amp, group_name, mice_colors, '.', SCATTER_CIRCLE_SIZE/2);
xlabel('Max amplitude of average ABR without stimulus (nV)')
ylabel(['Wave 1 amplitude @ ', num2str(THIS_INTENSITY), ' dB (nV)'])
title(['Wave 1 amplitude w/wo stimulus @ ', num2str(THIS_FREQ), ' Hz'])
set(gca,'FontSize', FONTSIZE) % Set text size of current axis

% Scatter plot of average ABR RMS with stimulus vs noise max amp
nexttile(t4)
gscatter(grouped_max_amp_noise_abr, grouped_rms_average_abr, group_name, mice_colors, '.', SCATTER_CIRCLE_SIZE/2);
xlabel('Max amplitude of average ABR without stimulus (nV)')
ylabel(['RMS of average ABR @ ', num2str(THIS_INTENSITY), ' dB (nV)'])
title(['RMS w stimulus vs peak-peak amp wo @ ', num2str(THIS_FREQ), ' Hz'])
set(gca,'FontSize', FONTSIZE) % Set text size of current axis

t4.TileSpacing = 'tight';
t4.Padding = 'tight';

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