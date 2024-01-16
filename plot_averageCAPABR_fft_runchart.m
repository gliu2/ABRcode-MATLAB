% function plot_averageCAPABR_fft_runchart(varargin)
% Plot amplitude of FFT of average CAP or ABR at specific
% noise frequency for all runs in selected folder.
%
% Default is for user to select file in dialog box when no input parameters
% are specified.
%
% 7/13/2023 George Liu
% Dependencies: import_ABRcsv.m, merge_singletraceABR_polarities, savefig_multiformat.m

%% Constants
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures
% SAMPLE_PERIOD = 40.96; % sample period of 40.96 us/sample in average trace
SAMPLE_PERIOD = 5.12; % sample period of 5.12 us/sample in single trace
SAMPLE_PERIOD_MS = SAMPLE_PERIOD/1000; % sample period in ms

NOISE_16KHZ = 8405; % Hz
NOISE_32KHZ = 7604; % Hz

NUM_CHANNELS = 2; % 1 = ABR, 2 = CAP
CHANNEL_KEY = {'ABR', 'CAP'};
CHANNEL_YLIM = {[-1000, 1000], [-1000, 1000]};

NUM_LEVELS = 10;
% N_FREQ = 2; % 16 kHz, 32 kHz

PLOT_FIGURES = 1; % Hard code logical for whether or not to plot figures 
SAVE_FIGURES = 1; % Hard code logical for whether or not to save figures 

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
dot_loc = strfind(filename, '.');

num_dashes = length(dash_loc);
group_ind = dash_loc(num_dashes - 3) + 1 : dash_loc(num_dashes - 2) - 1;
SGI_ind = dash_loc(num_dashes - 2) + 1 : dash_loc(num_dashes - 1) - 1;
channel_ind = dash_loc(num_dashes - 1) + 1 : dash_loc(num_dashes) - 1;
run_ind = dash_loc(num_dashes) + 1 : dot_loc - 1;

% Metadata of selected file
group = str2double(filename(group_ind));
SGI = str2double(filename(SGI_ind));
channel = str2double(filename(channel_ind)); % 1 = ABR, 2 = CAP
run = str2double(filename(run_ind));

% Obtain all single trial ABR/CAP files with same metadata (channel, SGI,
% run number)
filename_stem = filename(dash_loc(1):end);
listing = dir(fullfile(path, ['*', filename_stem]));

allnames = {listing.name}'; % column cell array of filenames

%% Find order of runs
NAME_STEM = 'testnoise_nomouse';
extractLabel = @(C, j, k) cellfun(@(c)  extractBetween(c, j, k), C);
extractLabelafter = @(C, k) cellfun(@(c)  extractAfter(c, k), C, 'UniformOutput',false);
allnames_labels = extractLabel(allnames, NAME_STEM, filename_stem);
allnames_labels2 = extractLabelafter(allnames_labels, '_');
allnames_ind = extractLabel(allnames, NAME_STEM, '_');
[~, run_ind] = sort(cellfun(@str2num, allnames_ind));

%%
num_data = length(run_ind);
runchart_fft = zeros(num_data, 1);
for i = 1:num_data
    this_filename = allnames{run_ind(i)};
    disp(['Working on file ', num2str(i), ' out of ', num2str(num_data), ': ', this_filename])
    [X_csv, A_csv, freq_csv] = import_ABRcsv(this_filename, path);

    % merge single trace pairs with alternating polarity to cancel cochlear
    % microphonic heterogeneity
    y = merge_singletraceABR_polarities(X_csv); % 1026 -> 513 single traces. Matrix of size (SAMPLES, m_traces/2)
    y_avg = mean(y, 2);

    num_samples = size(y, 1);
%         n_traces = size(y, 2);
    X = SAMPLE_PERIOD_MS * (1:num_samples); % time in ms

    % Calculate FFT value at noise frequency
    [f, P1] = get_fft_CAPABR(y_avg(1:num_samples/2));
    
    % Find value of Single-Sided Amplitude Spectrum (P1) nearest noise
    % frequency
    if freq_csv == 16000
        this_noise = NOISE_16KHZ;
        dist_noisefreq = abs(f - NOISE_16KHZ);
    elseif freq_csv == 32000
        this_noise = NOISE_32KHZ;
        dist_noisefreq = abs(f - NOISE_32KHZ);
    else
        disp('Warning: selected frequency does not have hardcoded noise frequency, skipping...')
        continue
    end
    [~, ind_freq] = min(dist_noisefreq);
    runchart_fft(i) =  P1(ind_freq);
end
        
%% Plot noise run chart 
fig = figure;
stem(1:num_data, runchart_fft, 'LineWidth', 3)

set(gca,'xtick',[1:num_data],'xticklabel', allnames_labels2(run_ind))
xtickangle(45)

% % Set text size of current axis
% set(gca,'FontSize',24)


%%
% xlabel('Amplitude (dB)', 'FontSize', 24)
ylabel(['|P1(f=', num2str(this_noise), ' Hz)| (nV)'], 'FontSize', 24)
title([CHANNEL_KEY{channel}, ' FFT run chart @ ', num2str(freq_csv), ' Hz, ', num2str(A_csv), ' dB'], 'FontSize', 24)

% Maximize figure window size
fig.WindowState = 'maximized';

if SAVE_FIGURES
    % Save figure
    %     disp('Saving figure')
    [~, save_file, ~] = fileparts(filename);
    savefig_multiformat(gcf, SAVE_PATH, [save_file, '_', num2str(A_csv), 'dB_', num2str(freq_csv/1000), 'kHz_', CHANNEL_KEY{channel}, '_fftrunchart2'])
end


disp('Done')

% end