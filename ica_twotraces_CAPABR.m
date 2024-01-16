% function ica_twotraces_CAPABR(varargin)
% Plot and determine ICA of average CAP and ABR at one stimulus level and frequency. Saves
% plots in multiple image formats.
%
% Default is user to select file in dialog box when no input parameters
% are specified.
%
% 7/19/2023 George Liu
% Dependencies: import_ABRcsv.m, merge_singletraceABR_polarities,
% savefig_multiformat.m, get_fft_CAPABR.m

opengl('save', 'software') % prevent crashing due to low-level graphics error using Sarah Office computer's graphics card

%% Constants
SAVE_PATH = 'd:\users\admin\Documents\George\Results_analyzed'; % path for saving figures
% SAMPLE_PERIOD = 40.96; % sample period of 40.96 us/sample in average trace
SAMPLE_PERIOD = 5.12; % sample period of 5.12 us/sample in single trace
SAMPLE_PERIOD_MS = SAMPLE_PERIOD/1000; % sample period in ms

CHANNEL_KEY = {'ABR', 'CAP'};
CHANNEL_COLOR = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980]}; % blue, red

N_ORDER = 700; % Half mainlobe width (-3 dB, approximation of cutoff frequency) of 226.484 Hz

FILTER_ON = 1; % Apply bandpass filter with Blackman window if no preacquisition filtering done
SAVE_FIGURES = 1; % Hard code logical for whether or not to save figures 

USE_ANALYZE_WINDOW = 1;
ANALYZE_WINDOW = [0, 3]; % ms; time window for analysis

%% Load average trace data
% if nargin == 2
%     path = varargin{1};
%     filename = varargin{2};
% elseif nargin == 0
%     [filename,path] = uigetfile('*.csv');
% else
%     disp('Warning: number of input arguments is not 0 or 2!')
%     [filename,path] = uigetfile('*.csv');
% end

[filename,path] = uigetfile('*.csv');
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

% Determine filename stem
filename_prestem = filename(1:dash_loc(3) - 1);

%% Load data
[X_csv, A_csv, freq_csv] = import_ABRcsv(filename, path);
if FILTER_ON 
    X_csv = filter_highpass_blackman(X_csv, N_ORDER); % Blackman high pass filter
end    

y = X_csv;
% y = merge_singletraceABR_polarities(X_csv); % 1026 -> 513 single traces. Matrix of size (SAMPLES, m_traces/2)
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

%% Load second average trace data (CAP if ABR chosen, or vice versa)
filename2 = filename;
if channel==1
    filename2(channel_ind) = '2';
    channel_plots = [1, 2]; % ABR then CAP
elseif channel==2
    filename2(channel_ind) = '1';
    channel_plots = [2, 1]; % CAP then ABR
end

% Load data
[X_csv2, A_csv2, freq_csv2] = import_ABRcsv(filename2, path);
if FILTER_ON 
    X_csv2 = filter_highpass_blackman(X_csv2, N_ORDER); % Blackman high pass filter
end
    
y2 = X_csv2;
% y = merge_singletraceABR_polarities(X_csv); % 1026 -> 513 single traces. Matrix of size (SAMPLES, m_traces/2)
y_avg2 = mean(y2, 2);

%% Plot overlay of two average traces (ABR and CAP)
% Calculate ylims for plotting
ylim_max = [];
ymax_data = max(y_avg, [], "all");
ymin_data = min(y_avg, [], "all");
ymax_data2 = max(y_avg2, [], "all");
ymin_data2 = min(y_avg2, [], "all");
ylim_max(2) = max([ymax_data, ymax_data2]);
ylim_max(1) = min([ymin_data, ymin_data2]);

fig = figure;
plot(X, y_avg, 'LineWidth', 2.5, 'Color', CHANNEL_COLOR{channel_plots(1)})
hold on
plot(X, y_avg2, 'LineWidth', 2.5, 'Color', CHANNEL_COLOR{channel_plots(2)})
hold off
ylim(ylim_max)
lg = legend(CHANNEL_KEY(channel_plots), 'Location', 'northeast');

% Set text size of current axis
set(gca,'FontSize',24)

% Show frequency as title above top-most tile of each column
title(['ABR and CAP @ ', num2str(freq_csv), ' Hz, ', num2str(A_csv), ' dB']) 
xlabel('Time (ms)', 'FontSize', 24)
ylabel('nV', 'FontSize', 24)

% Maximize figure window size
fig.WindowState = 'maximized';

if SAVE_FIGURES
    save_filename = [filename_prestem, '_', num2str(A_csv), 'dB_', num2str(freq_csv/1000), 'kHz_overlayCAPABR'];
    savefig_multiformat(gcf, SAVE_PATH, save_filename)
end

%% Reduce window for analysis
if USE_ANALYZE_WINDOW 
    n_coords = length(ANALYZE_WINDOW);
    ind_window = zeros(n_coords, 1);
    for i = 1:n_coords
        dist = abs(X - ANALYZE_WINDOW(i));
        [~, ind_window(i)] = min(dist);
    end
    
    X = X(ind_window(1):ind_window(2));
    y_avg = y_avg(ind_window(1):ind_window(2));
    y_avg2 = y_avg2(ind_window(1):ind_window(2));
end

%% Zero mean and whiten data prior to ICA analysis
% y_mean = mean([y_avg, y_avg2], 2);
% y_std = std([y_avg, y_avg2], 0, 2);
% 
% y_avg = (y_avg - y_mean) ./ y_std;
% y_avg2 = (y_avg2 - y_mean) ./ y_std;

y_avg = (y_avg - mean(y_avg)) / rms(y_avg);
y_avg2 = (y_avg2 - mean(y_avg2)) / rms(y_avg2);
    
%% Calculate ICA of average CAP and ABR traces, at single frequency and stimulus 
rng default
q = 2;
Mdl = rica([y_avg'; y_avg2'], q); % input is matrix of row vectors

%% Plot independent components
% Apply ICA
% data_ICA = transform(Mdl, y');
comps_ICA = Mdl.TransformWeights;

NUM_COLS = 1;
fig3 = figure;
t2 = tiledlayout(q/NUM_COLS , NUM_COLS, 'TileIndexing', 'rowmajor');
for i = 1:q
    nexttile(t2)
    plot(X, comps_ICA(:, i), 'LineWidth', 3, 'Color', [0, 0, 0])
    title(strcat("Component ", string(i)))
    
    hold on
    plot(X, y_avg/norm(y_avg), 'LineWidth', 2.5, 'Color', CHANNEL_COLOR{channel_plots(1)}, 'LineStyle', '--')
    plot(X, y_avg2/norm(y_avg2), 'LineWidth', 2.5, 'Color', CHANNEL_COLOR{channel_plots(2)}, 'LineStyle', '--')
    hold off
    
    ylim([-0.1, 0.1])
    set(gca,'FontSize',24)
end
lg2 = legend({'Comp', CHANNEL_KEY{channel_plots(1)}, CHANNEL_KEY{channel_plots(2)}}, 'Orientation', 'Vertical');
lg2.Layout.Tile = 'East'; % <-- Legend placement with tiled layout

% Set text size of current axis
set(gca,'FontSize',24)

title(t2, ['ICA @ ', num2str(freq_csv), ' Hz, ', num2str(A_csv), ' dB'], 'FontSize', 24)
xlabel(t2, 'Time (ms)', 'FontSize', 24)
ylabel(t2, 'Normalized amplitude', 'FontSize', 24)

t2.TileSpacing = 'tight';
t2.Padding = 'tight';

fig3.WindowState = 'maximized'; % Maximize figure window size

if SAVE_FIGURES
    if USE_ANALYZE_WINDOW
        save_filename2 = [filename_prestem, '_', num2str(A_csv), 'dB_', num2str(freq_csv/1000), 'kHz_CAPABR_ICA_', num2str(q), 'comps_window'];
    else
        save_filename2 = [filename_prestem, '_', num2str(A_csv), 'dB_', num2str(freq_csv/1000), 'kHz_CAPABR_ICA_', num2str(q), 'comps'];
    end
    savefig_multiformat(gcf, SAVE_PATH, save_filename2)
end
disp('Done')

% end