function plot_average_CAPmicrophonic(varargin)
% Plot average and microphonic of single traces file of CAP or ABR at one stimulus level and frequency. Saves
% plots in multiple image formats.
%
% Default is user to select file in dialog box when no input parameters
% are specified.
%
% 7/10/2023 George Liu
% Dependencies: import_ABRcsv.m, merge_singletraceABR_polarities, savefig_multiformat.m

%% Constants
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures
% SAMPLE_PERIOD = 40.96; % sample period of 40.96 us/sample in average trace
SAMPLE_PERIOD = 5.12; % sample period of 5.12 us/sample in single trace
SAMPLE_PERIOD_MS = SAMPLE_PERIOD/1000; % sample period in ms

CHANNEL_KEY = {'ABR', 'CAP'};
YLIM_DEFAULT = [-1000, 1000];

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
    
% merge single trace pairs with alternating polarity to cancel cochlear
% microphonic heterogeneity
y = merge_singletraceABR_polarities(X_csv); % 1026 -> 513 single traces. Matrix of size (SAMPLES, m_traces/2)
y_avg = mean(y, 2);

num_samples = size(y, 1);
X = SAMPLE_PERIOD_MS * (1:num_samples); % time in ms

% % TESTING: use single frequency trace as simulated data for testing fft
% F_test = 16; % kHz
% y_avg = sin(2*pi*F_test*X);

% Make sure bounds of plot are within y limits
ylim_max = YLIM_DEFAULT;
ymax_data = max(y_avg);
ymin_data = min(y_avg);
if ymax_data > ylim_max(2)
    ylim_max(2) = ymax_data;
end
if ymin_data < ylim_max(1)
    ylim_max(1) = ymin_data;
end

fig = figure;
t = tiledlayout(1, 3);
nexttile(t)
plot(X, y_avg, 'LineWidth', 3)
ylim(ylim_max)

ylabel([num2str(A_csv), ' dB (nV)'], 'FontSize', 24)
% rotate y label to make horizontal
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right', 'FontSize', 24)

% Adjust plot appearance
set(gca,'box','off')

% Set text size of current axis
set(gca,'FontSize',24)

% Show frequency as title above top-most tile of each column
title([CHANNEL_KEY{channel}, ' @ ', num2str(freq_csv), ' Hz']) 
xlabel('Time (ms)', 'FontSize', 24)


%% Plot cochlear microphonic
nexttile(t)
y_microphonic = subtract_singletraceABR_polarities(X_csv);
y_avg2 = mean(y_microphonic, 2);
plot(X, y_avg2, 'LineWidth', 3)
title(['Microphonic @ ', num2str(freq_csv), ' Hz']) 
xlabel('Time (ms)', 'FontSize', 24)

ylabel("nV")

% Adjust plot appearance
% Remove extraneous axis tick marks and x-axis from all but bottom
% tile
set(gca,'box','off')

% Set text size of current axis
set(gca,'FontSize',24)

%% Plot FFT of microphonic trace
y_fft = fft(y_avg2);
Fs = 1000/SAMPLE_PERIOD_MS; % Sampling frequency 
L = length(y_avg); % sample length 
P2 = abs(y_fft/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;


nexttile(t)
plot(f,P1,'LineWidth', 3) 
xlim([0, 33000])

% FFT plot labels
title("1-Sided FFT microphonic")
xlabel("f (Hz)")
ylabel("|P1(f)|")

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
    savefig_multiformat(gcf, SAVE_PATH, [save_file, '_', num2str(A_csv), 'dB_', num2str(freq_csv/1000), 'kHz_', CHANNEL_KEY{channel}, '_microphonic'])
end

disp('Done')

end