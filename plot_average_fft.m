function plot_average_fft(varargin)
% Plot FFT of selected single traces file of CAP or ABR at one stimulus level and frequency. Saves
% plots in multiple image formats.
%
% Default is user to select file in dialog box when no input parameters
% are specified.
%
% 8/20/2023 George Liu
% Dependencies: notch_filter.m, get_fft_CAPABR.m, import_ABRcsv.m, merge_singletraceABR_polarities, savefig_multiformat.m

%% Constants
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures
% SAMPLE_PERIOD = 40.96; % sample period of 40.96 us/sample in average trace
SAMPLE_PERIOD = 5.12; % sample period of 5.12 us/sample in single trace
SAMPLE_PERIOD_MS = SAMPLE_PERIOD/1000; % sample period in ms

CHANNEL_KEY = {'ABR', 'CAP'};
YLIM_DEFAULT = [-1000, 1000];
XLIM_MAX = 10000;

N_ORDER = 700; % Half mainlobe width (-3 dB, approximation of cutoff frequency) of 226.484 Hz
FILTER_HIGHPASS_ON = 1; % Apply high pass filter with Blackman window of N_ORDER length. Use if no preacquisition high pass filtering was done.

SAVE_FIGURES = 1; % Hard code logical for whether or not to save figures 

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

% Analyze metadata in single trial ABR/CAP filename
% Filename example: 20230628_tmie2_ABRCAP_pre-0-55-2-1.csv
% group 0
% SGI 55
% channel 2
% run (?) 1
dash_loc = strfind(filename, '-');

num_dashes = length(dash_loc);
channel_ind = dash_loc(num_dashes - 1) + 1 : dash_loc(num_dashes) - 1;
channel = str2double(filename(channel_ind));

% Load data
[X_csv, A_csv, freq_csv] = import_ABRcsv(filename, path);

% Filter data
if FILTER_HIGHPASS_ON 
    X_csv = filter_highpass_blackman(X_csv, N_ORDER); % Blackman high pass filter  
    
    % Apply notch filter to remove 1901.1 Hz noise
    X_csv = notch_filter(X_csv);
    
end

% merge single trace pairs with alternating polarity to cancel cochlear
% microphonic heterogeneity
y = merge_singletraceABR_polarities(X_csv); % 1026 -> 513 single traces. Matrix of size (SAMPLES, m_traces/2)
y_avg = mean(y, 2);

num_samples = size(y, 1);
X = SAMPLE_PERIOD_MS * (1:num_samples); % time in ms

% % Make sure bounds of plot are within y limits
% ylim_max = YLIM_DEFAULT;
% ymax_data = max(y_avg);
% ymin_data = min(y_avg);
% if ymax_data > ylim_max(2)
%     ylim_max(2) = ymax_data;
% end
% if ymin_data < ylim_max(1)
%     ylim_max(1) = ymin_data;
% end
% 
% fig = figure;
% t = tiledlayout(1, 2);
% nexttile(t)
% plot(X, y_avg, 'LineWidth', 3)
% ylim(ylim_max)
% 
% ylabel([num2str(A_csv), ' dB (nV)'], 'FontSize', 24)
% % rotate y label to make horizontal
% hYLabel = get(gca,'YLabel');
% set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right', 'FontSize', 24)
% 
% % Adjust plot appearance
% set(gca,'box','off')
% 
% % Set text size of current axis
% set(gca,'FontSize',24)
% 
% % Show frequency as title above top-most tile of each column
% title([CHANNEL_KEY{channel}, ' @ ', num2str(freq_csv), ' Hz']) 
% xlabel('Time (ms)', 'FontSize', 24)
% 
% % Obtain FFT of stack plots
% [f, P1] = get_fft_CAPABR(y_avg);
% 
% nexttile(t)
% plot(f,P1,'LineWidth', 3) 
% xlim([0, XLIM_MAX])
% 
% % FFT plot labels
% title("Single-Sided Amplitude Spectrum of X(t)")
% xlabel("f (Hz)")
% ylabel("Amplitude (nV)")
% 
% % Adjust plot appearance
% % Remove extraneous axis tick marks and x-axis from all but bottom
% % tile
% set(gca,'box','off')
% 
% % Set text size of current axis
% set(gca,'FontSize',24)
% 
% title(t, filename, 'FontSize', 24, 'Interpreter', 'none')
% 
% % Maximize figure window size
% fig.WindowState = 'maximized';
% 
% if SAVE_FIGURES
%     % Save figure
%     %     disp('Saving figure')
%     [~, save_file, ~] = fileparts(filename);
%     savefig_multiformat(gcf, SAVE_PATH, [save_file, '_', num2str(A_csv), 'dB_', num2str(freq_csv/1000), 'kHz_', CHANNEL_KEY{channel}, '_fft'])
% end

%% Plot CAP FFT, odd traces FFT, and even traces FFT
y_cm = merge_singletraceABR_polarities_CM(X_csv); % 1026 -> 513 single traces. Matrix of size (SAMPLES, m_traces/2)
y_avg_cm = mean(y_cm, 2);
y_avg_odd = (y_avg + y_avg_cm)/2 ;
y_avg_even = (y_avg - y_avg_cm)/2 ;

data = {y_avg, y_avg_cm, y_avg_odd, y_avg_even};
n_data = length(data);
data_labels = {'', 'CM', 'odd', 'even'};

for i = 1:n_data
% for i = 1:1
    this_data = data{i};
    
    fig = figure;
    t = tiledlayout(1, 2);
    nexttile(t)
    plot(X, this_data, 'LineWidth', 3)
    ylim([min(this_data), max(this_data)])

    ylabel('Amplitude (nV)', 'FontSize', 24)
    
    % Adjust plot appearance
    set(gca,'box','off')

    % Set text size of current axis
    set(gca,'FontSize',24)

    % Show frequency as title above top-most tile of each column
    title([CHANNEL_KEY{channel}, ' ', data_labels{i}, ' @ ', num2str(freq_csv), ' Hz, ', num2str(A_csv), ' dB']) 
    xlabel('Time (ms)', 'FontSize', 24)

    % Obtain FFT of stack plots
    [f, P1, phi] = get_fft_CAPABR(this_data);

    nexttile(t)
    yyaxis left
    plot(f,P1,'LineWidth', 3) 
    ylabel("Single-sided amplitude (nV)")
    ylim([0, 1500])
    xlim([0, XLIM_MAX])
    
    yyaxis right
    plot(f, unwrap(phi), 'LineWidth', 2, 'LineStyle', ':')
%     ylim([-pi, pi])
    ylabel("Phase (rad)")
    
    xline(8414, '-', 'Aliased 16 kHz signal', 'FontSize', 16)
    xline(7586, '-', 'Aliased 32 kHz signal', 'FontSize', 16)
    xline(8000, '-', '8 kHz', 'FontSize', 16)
    xline(1901.1, '-', '1901.1 Hz', 'FontSize', 16)

    % FFT plot labels
    title("Fourier Transform")
    xlabel("f (Hz)")
    

    % Adjust plot appearance
    % Remove extraneous axis tick marks and x-axis from all but bottom
    % tile
    set(gca,'box','off')

    % Set text size of current axis
    set(gca,'FontSize',24)

    title(t, filename, 'FontSize', 24, 'Interpreter', 'none')

    % Maximize figure window size
    fig.WindowState = 'maximized';

    if SAVE_FIGURES
        % Save figure
        %     disp('Saving figure')
        [~, save_file, ~] = fileparts(filename);
        savefig_multiformat(gcf, SAVE_PATH, [save_file, '_', num2str(A_csv), 'dB_', num2str(freq_csv/1000), 'kHz_', CHANNEL_KEY{channel}, '_', data_labels{i}, '_fft'])
    end

disp('Done')

end