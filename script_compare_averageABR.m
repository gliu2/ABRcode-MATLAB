% script_compare_averageABR
%
% 7/23/21 George Liu
% Dependencies: parse_mouse_name.m, get_mousegroup_fromkey.m

%% Constants
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures
PATH_MOUSEKEY = 'd:\users\admin\Documents\George\mice_key.xlsx';
% PATH_NOISE_RESULTS = 'd:\users\admin\Documents\George\aggregate_results\noise_exposure_data.xlsx';
PATH_NOISE_RESULTS = 'd:\users\admin\Documents\George\aggregate_results\control_data.xlsx';

%% Load data
T = readtable(PATH_NOISE_RESULTS, 'Range','A3');

% cols = {'ABRThreshold', 'ABRThreshold_1', 'ABRThreshold_2', 'ABRThreshold_3', 'ABRThreshold_4', 'ABRThreshold_5'};
% cols = {'DPOAEThreshold', 'DPOAEThreshold_1', 'DPOAEThreshold_2', 'DPOAEThreshold_3', 'DPOAEThreshold_4', 'DPOAEThreshold_5'};
cols = {'Amplitude', 'Amplitude_1', 'Amplitude_2', 'Amplitude_3', 'Amplitude_4', 'Amplitude_5'};


%% Plot threshold shift
pre_ABRthresh_8khz = T.(cols{1});
pre_ABRthresh_16khz = T.(cols{2});
pre_ABRthresh_32khz = T.(cols{3});
pre_ABRthresh = horzcat(pre_ABRthresh_8khz, pre_ABRthresh_16khz, pre_ABRthresh_32khz);
post_ABRthresh_8khz = T.(cols{4});
post_ABRthresh_16khz = T.(cols{5});
post_ABRthresh_32khz = T.(cols{6});
post_ABRthresh = horzcat(post_ABRthresh_8khz, post_ABRthresh_16khz, post_ABRthresh_32khz);

freq = [8, 16, 32];
% jitter = rand(size(freq))*2;
figure
errorbar(freq, nanmean(pre_ABRthresh, 1), nanstd(pre_ABRthresh, 1), 'b')
hold on
errorbar(freq, nanmean(post_ABRthresh, 1), nanstd(post_ABRthresh, 1), 'r')
hold off

% ylim([0, 90])
xlabel('Frequency (kHz)')
% ylabel('Threshold (dB)')
ylabel('Amplitude (nV)')
% title('Control')
% title('Noise exposure')
% title('Wave 1 amplitude, noise')
title('Wave 1 amplitude, control')
legend('Pre-noise', '24h post')