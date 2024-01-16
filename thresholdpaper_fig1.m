% function thresholdpaper_fig1(varargin)
% Plot of example average CAP waveforms, single trial inner product
% histograms, D' calculation, and threshold algorithm based on D'.
%
% Last edit: 8/29/23 George Liu
%
% Dependencies: get_thresh_CAPABR_D.m, steshade.m, stdshade.m

PATH_ROOT = '\\Ricci-abr\d\George\CAP_ABR\';
PATH_ALLCONTROLS = {...
%     Control no tube - pre
    '20230823\20230823_control1_pre-0-24-1-1.csv'; ...
    '20230823\20230823_control2_pre-0-24-1-1.csv'; ... % no baseline hearing at 32 kHz
    '20230823\20230823_control3_pre-0-24-1-1.csv'; ...
    '20230823\20230823_control4_pre-0-24-1-1.csv'; ...
    '20230824\20230824_control1_pre-0-24-1-1.csv'; ...
    %     Control no tube - pre2
    '20230823\20230823_control1_pre2-0-24-1-1.csv'; ...
    '20230823\20230823_control2_pre2-0-24-1-1.csv'; ... % no baseline hearing at 32 kHz
    '20230823\20230823_control3_pre2-0-24-1-1.csv'; ...
    '20230823\20230823_control4_pre2-0-24-1-1.csv'; ...
    '20230824\20230824_control1_pre2-0-24-1-1.csv'; ...
    %     Control no tube - pre3
    '20230823\20230823_control1_pre3-0-24-1-1.csv'; ...
    '20230823\20230823_control2_pre3-0-24-1-1.csv'; ... % no baseline hearing at 32 kHz
    '20230823\20230823_control3_pre3-0-24-1-1.csv'; ...
    '20230823\20230823_control4_pre3-0-24-1-1.csv'; ...
    '20230824\20230824_control1_pre3-0-24-1-1.csv'; ...
    %     Control no tube - pre4
    '20230823\20230823_control1_pre4-0-24-1-1.csv'; ...
%     '20230823\20230823_control2_pre4-0-24-1-1.csv'; ... % doesn't exist
    '20230823\20230823_control3_pre4-0-24-1-1.csv'; ...
    '20230823\20230823_control4_pre4-0-24-1-1.csv'; ...
    '20230824\20230824_control1_pre4-0-24-1-1.csv' ...
};
% PATH_ALLCONTROLS = {...
% %     Control tube - pre
%     '20230725\20230725_control1_pre-0-24-1-1.csv'; ...
%     '20230726\20230726_control1_pre-0-24-1-1.csv'; ... % no baseline hearing at 32 kHz
%     '20230726\20230726_control2_pre-0-24-1-1.csv'; ...
%     '20230828\20230828_tube2_pre-0-24-1-1.csv'; ...
%     % Control AP - pre
%     '20230727\20230727_control1_pre-0-24-1-1.csv'; ...
%     '20230727\20230727_control2_pre-0-24-1-1.csv'; ...
%     '20230727\20230727_control3_pre-0-24-1-1.csv'; ... % no baseline hearing at 32 kHz
%     '20230813\20230813_control2_pre-0-24-1-1.csv'; ...
%     '20230903\20230903_control1_pre-0-24-1-1.csv'; ... % Exclude because waves change shape significantly
%     '20230903\20230903_control2_pre-0-24-1-1.csv'; ...
%     '20230903\20230903_control3_pre-0-24-1-1.csv'; ... % Exclude because waves change shape significantly
%     % CNQX_FILES 
%     '20230801\20230801_cnqx1_pre-0-24-1-1.csv'; ... % no baseline hearing at 32 kHz
%     '20230801\20230801_cnqx2_pre-0-24-1-1.csv'; ...
%     '20230802\20230802_cnqx1_pre-0-24-1-1.csv'; ...
%     '20230828\20230828_cnqx1_pre-0-24-1-1.csv'; ...
%     '20230828\20230828_cnqx3_pre-0-24-1-1.csv'; ...
%     '20230829\20230829_cnqx1_pre-0-24-1-1.csv'; ...
%     '20230829\20230829_cnqx2_pre-0-24-1-1.csv'; ...
%     '20230831\20230831_cnqx1_pre-0-24-1-1.csv'; ...
%     '20230831\20230831_cnqx2_pre-0-24-1-1.csv'; ...
%     '20230831\20230831_cnqx3_pre-0-24-1-1.csv'; ...
% %     BAD_AP_FILES 
%     '20230721\20230721_control1_pre-0-24-1-1.csv'; ...
%     '20230722\20230722_control1_pre-0-24-1-1.csv'; ...
%     '20230723\20230723_control1_pre-0-24-1-1.csv'; ... 
%     '20230724\20230724_control1_pre-0-24-1-1.csv'; ... 
%     '20230724\20230724_control2_pre-0-24-1-1.csv' ... 
% };
% PATH_ALLCONTROLS = {...
%     '20230724\20230724_control2_pre-0-24-1-1.csv' ... 
% };

% PATH_EXAMPLE_CAP = '\\Ricci-abr\d\George\CAP_ABR\20230823\20230823_control3_pre-0-24-1-1.csv';
PATH_EXAMPLE_CAP = '\\Ricci-abr\d\George\CAP_ABR\20230813\20230813_control2_pre-0-24-1-1.csv';
% PATH_EXAMPLE_CAP = '\\Ricci-abr\d\George\CAP_ABR\20230726\20230726_control2_pre-0-24-1-1.csv';
% PATH_EXAMPLE_CAP = '\\Ricci-abr\d\George\CAP_ABR\20230801\20230801_cnqx1_pre-0-24-1-1.csv';
% SAVE_PATH = 'd:\users\admin\Documents\George\Results_analyzed'; % path for saving figures
SAVE_PATH = 'd:\users\admin\Documents\George\Results_analyzed_10-23-23'; % for no tube control data analysis

SAMPLE_PERIOD = 5.12; % sample period of 5.12 us/sample in single trace
SAMPLE_PERIOD_MS = SAMPLE_PERIOD/1000; % sample period in ms
NV_TO_UV_CONVERSION = 10^-3;

myColors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980]}; % blue and orange
FONT_SIZE = 24;
LINE_WIDTH = 2;

THIS_FREQUENCY = [8000, 16000, 32000];
% THIS_FREQUENCY = 16000; % Frequency stack to plot
% THIS_FREQUENCY = 32000; % Frequency stack to plot
N_COLUMNS = 3;

NUM_CHANNELS = 2; 
THESE_CHANNELS = 1:2; % 1 = ABR, 2 = CAP
% THESE_CHANNELS = 1; % 1 - ABR only; 2 - CAP only
CHANNEL_KEY = {'ABR', 'ECochG'};
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

PLOT_WAVE1 = 0; % Plot overlay circles on peak and trough of wave 1
SAVE_FIGURES = 1; % Hard code logical for whether or not to save figures 

%% Begin analysis
num_experiments = length(PATH_ALLCONTROLS);
for aa = 1:num_experiments
    % Close existing windows to avoid crashing
    close all
    
    % Get this experiment's filename(s)
    this_filepath = fullfile(PATH_ROOT, PATH_ALLCONTROLS{aa});

%     % [filename,path] = uigetfile('*.csv'); % Uncomment if ui select file
%     [path, filename, ext] = fileparts(PATH_EXAMPLE_CAP);
%     filename = [filename, ext];
%     disp(['Opening ', PATH_EXAMPLE_CAP])
    [path, filename, ext] = fileparts(this_filepath);
    filename = [filename, ext];
    disp(['  Opening ', num2str(aa), ' out of ', num2str(num_experiments), ': ', this_filepath])

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

    for i = THESE_CHANNELS
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

        % Collect all average ABR traces from single trial ABR recordings
        y_avg_cache = cell(num_data, 1);
        y_avg_cm_cache = cell(num_data, 1); % microphonic
        A_csv_cache = zeros(num_data, 1);
        freq_csv_cache = zeros(num_data, 1);
        X_csv_cache = cell(num_data, 1);

        for j = 1:num_data
    %     for j = 11:20
            this_filename = data_filenames_ordered{j};
            disp(['Working on file ', num2str(j), ' out of ', num2str(num_data), ': ', this_filename])
            [X_csv, A_csv, freq_csv] = import_ABRcsv(this_filename, path);

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

            % Convert nV to uV
            X_csv = X_csv * NV_TO_UV_CONVERSION;

            % merge single trace pairs with alternating polarity to cancel cochlear
            % microphonic heterogeneity
            y = merge_singletraceABR_polarities(X_csv); % 1026 -> 513 single traces. Matrix of size (SAMPLES, m_traces/2)
            y_avg = mean(y, 2);

            num_samples = size(y, 1);
        %    n_traces = size(y, 2);
            X = SAMPLE_PERIOD_MS * (1:num_samples); % time in ms

            % Cache variables
            y_avg_cache{j, 1} = y_avg'; % make row vector
            A_csv_cache(j, 1) = A_csv;
            freq_csv_cache(j, 1) = freq_csv;
            X_csv_cache{j, 1} = X_csv;

        end

        %% Reshape cache of single trial data
        X_csv_cache = reshape(X_csv_cache, [NUM_LEVELS, n_freq]);
        freq_csv_cache = reshape(freq_csv_cache, [NUM_LEVELS, n_freq]);
        A_csv_cache = reshape(A_csv_cache, [NUM_LEVELS, n_freq]);

        %% Iterate over all frequencies
        for ii = 1:length(THIS_FREQUENCY)
            this_freq = THIS_FREQUENCY(ii);

            % Calculate threshold using single trial inner products
            ind_freq = find(freq_csv_cache(1, :)==this_freq , 1);
            [A_descending, I] = sort(A_csv_cache(:, ind_freq), 'descend');
            cutoff = 1;
            [threshold_D, metric_D, dist_innerprod] = get_thresh_CAPABR_D(X_csv_cache(:, ind_freq), A_descending, cutoff, false);

            %% Plot average traces in tiled plot
            fig = figure;
            t = tiledlayout(NUM_LEVELS, N_COLUMNS, 'TileIndexing', 'columnmajor');

            % Make sure bounds of plot are within y limits
            ymax_data = max(cell2mat(y_avg_cache), [], "all");
            ymin_data = min(cell2mat(y_avg_cache), [], "all");
            ylim_max = [ymin_data, ymax_data];

            for j = 1:num_data
                y_avg = y_avg_cache{j};
                A_csv = A_csv_cache(j);
                freq_csv = freq_csv_cache(j);

                % Only plot chosen frequency
                if freq_csv ~= this_freq  
                    continue
                end

                nexttile(t)
            %     plot(X, y_avg, 'LineWidth', 3)
                stdshade(X_csv_cache{j}', 0.2, myColors{1}, X, [], 3);
            %     steshade(X_csv_cache{j}', 0.1, myColors{1}, X, [], 3);
%                 ylim(ylim_max)
                ylim([-2.5, 2.5])
                yticks([-2, 0, 2])
%                 yticks([0, 2])

                % show ylabels for first column only
                ylabel([num2str(A_csv), ' dB (\muV)'], 'FontSize', FONT_SIZE)
                % rotate y label to make horizontal
                hYLabel = get(gca,'YLabel');
                set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right', 'FontSize', FONT_SIZE)

                % Adjust plot appearance
                % Remove extraneous axis tick marks and x-axis from all but bottom
                % tile
                set(gca,'box','off')
                if mod(j, NUM_LEVELS) ~= 0 % 
                    set(gca,'xtick',[])
                    set(gca,'xticklabel',[])
                    set(gca,'XColor','none')
                end

                % Set text size of current axis
                set(gca,'FontSize', FONT_SIZE)

                % Show title of column one
                if mod(j, NUM_LEVELS)==1
                    title(['Mouse ', CHANNEL_KEY{i}]) 
                end

                xlabel('Time (ms)', 'FontSize', FONT_SIZE)
            %     xlim([0, 12.5])
            end

            % Plot distributions of inner products
            max_innerprod = max(dist_innerprod, [], "all");
            min_innerprod = min(dist_innerprod, [], "all");
            innerprod_xlim = [min_innerprod, max_innerprod];
            max_bin_val = zeros(NUM_LEVELS, 1);
            axesHandle = zeros(NUM_LEVELS, 1);

            for j = 1:NUM_LEVELS
                axesHandle(j) = nexttile(t);
                % Plot histogram of single-trace inner products
                signal_mean = mean(dist_innerprod(j, :));
            %     h = histogram(dist_innerprod(j, :), 'BinMethod', 'fd', 'Normalization', 'probability',  'LineWidth', 0.5, 'FaceAlpha', 0.5, 'edgecolor', 'none');
                hold on    
                if j < NUM_LEVELS
                    h = histogram(dist_innerprod(j, :), 'BinMethod', 'fd', 'Normalization', 'probability',  'DisplayStyle','stairs', 'LineWidth', LINE_WIDTH, 'edgecolor', myColors{1});
                    xline(signal_mean, 'LineWidth', LINE_WIDTH, 'Color', myColors{1})
                else
                    h = histogram(dist_innerprod(j, :), 'BinMethod', 'fd', 'Normalization', 'probability',  'DisplayStyle','stairs', 'LineWidth', LINE_WIDTH, 'edgecolor', myColors{2});
                    xline(signal_mean, 'LineWidth', LINE_WIDTH, 'Color', myColors{2})
                end

                xlabel('Inner product (\muV^2)')

                % plot mean as blue line
            %     line([signal_mean signal_mean], ylim, 'LineWidth', 1, 'Color', 'b'); % vertical line for cutoff
                hold off

                if j==1
                    title('Histogram') 
                end

                % show ylabels for first column only
                if j==5
                    ylabel('Density', 'FontSize', FONT_SIZE)
            %     % rotate y label to make horizontal
            %     hYLabel = get(gca,'YLabel');
            %     set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right', 'FontSize', 24)
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

                % Set text size of current axis
                set(gca,'FontSize', FONT_SIZE)

                % Same x axes for all histograms
                xlim(innerprod_xlim)
                max_bin_val(j) = max(h.Values, [], "all");

                % Set y axis to log
                set(gca,'YScale','log')
            end
            % Set y axes same
            set(axesHandle, 'YLim', [0, max(max_bin_val)]);

            % Plot D'
            nexttile([NUM_LEVELS, 1])
            plot(metric_D, A_descending, 'LineWidth', 3)
            xlabel('Cohen''s d', 'FontSize', FONT_SIZE)
            my_xline = xline(cutoff, '--', ['Cutoff @ d = ', num2str(cutoff)], 'LineWidth', LINE_WIDTH);
            % title(['Threshold at ', num2str(round(threshold_D, 1)), ' dB'], 'FontSize', FONT_SIZE)
            title(['Threshold'], 'FontSize', FONT_SIZE)
            ylabel('Stimulus level (dB SPL)', 'FontSize', FONT_SIZE)
            set(my_xline, 'FontSize', FONT_SIZE)
            set(gca,'TickDir','out');
            set(gca,'box','off')
            % Set text size of current axis
            set(gca,'FontSize', FONT_SIZE)
            ylim([A_descending(end)-5, A_descending(1)+5]);
            xlim_bounds = [min([metric_D; 0]), 1.05*max(metric_D)];
            xlim(xlim_bounds);

            % % Set x axis to log
            % set(gca,'XScale','log')

            % Draw horizontal line at threshold amplitude 
            if ~isnan(threshold_D)
                hold on
                my_yline = yline(threshold_D, '--', ['Threshold = ', num2str(round(threshold_D, 1)), ' dB SPL'], 'LineWidth', LINE_WIDTH, 'Color', 'k'); % horizontal line at exact threshold
                set(my_yline, 'FontSize', FONT_SIZE)
                hold off
            else
                disp('Warning: No ABR threshold detected!')
            end

            % Format figure
            % xlabel(t, 'Time (ms)', 'FontSize', FONT_SIZE)
            ylabel(t, 'Stimulus level (dB SPL)', 'FontSize', FONT_SIZE)

            % t.TileSpacing = 'none';
            t.TileSpacing = 'tight';
            t.Padding = 'tight';

            % Maximize figure window size
            fig.WindowState = 'maximized';    

            % Create title
            figure_title = [filename_stem, ' - ', CHANNEL_KEY{i}, ', ', num2str(this_freq/1000), ' kHz'];
            title(t, figure_title, 'FontSize', FONT_SIZE, 'Interpreter','none')

            %% Save figure
            if SAVE_FIGURES
                [~, save_file, ~] = fileparts(filename_stem);
                savefig_multiformat(gcf, SAVE_PATH, [save_file, '_', CHANNEL_KEY{i}, '_freq_', num2str(this_freq), '_thresholdpaper_Fig1'])
            end
        end
    end
end

disp('Done')