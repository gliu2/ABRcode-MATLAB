function view_microphonic_CAPABR(varargin)
% Plot microphinic of average CAP and ABR data 
%
% 8/18/2023 George Liu
% Dependencies: tileplot_CAPABR.m

%% Constants
LOAD_PATH = 'd:\users\admin\Documents\George\Results';
SAVE_PATH = 'd:\users\admin\Documents\George\Results_analyzed'; % path for saving figures

SAMPLE_PERIOD = 5.12; % sample period of 5.12 us/sample in single trace
SAMPLE_PERIOD_MS = SAMPLE_PERIOD/1000; % sample period in ms

NUM_CHANNELS = 2; % 1 = ABR, 2 = CAP
CHANNEL_KEY = {'ABR', 'CAP'};

%% Load data
disp('Select MAT file in Results folder: ')
[filename,path] = uigetfile('*.mat', 'Select MAT file');
disp(['Opening ', fullfile(path, filename)])

% Determine name of pre injection MAT file
underscore_loc = strfind(filename, '_');

n_underscore = length(underscore_loc);
substrings = cell(n_underscore, 1);
substrings{1} = filename(1:underscore_loc(1)-1);
for i=1:n_underscore-1
    substrings{i+1} = extractBetween(filename, underscore_loc(i), ...
        underscore_loc(i+1), 'Boundaries', 'exclusive');
end

% date = substrings{1};
% label = substrings{2}{1};
% timepoint = substrings{3}{1};
% % ABRCAP = substring{4}{1};
% 
% filename_prestem = [date, '_', label, '_'];
my_vars = {'y_avg_cache_descending', 'y_avg_cm_cache_descending', ...
    'all_freqs', 'A_descending', 'threshold', 'metric'};

i = 2; % CAP only
% filename_poststem = ['_', CHANNEL_KEY{i}, '_results.mat'];
% %     filename_pre = [filename_prestem, 'pre', filename_poststem];
% filename_post = [filename_prestem, timepoint, filename_poststem];

% Load selected and baseline ABR / CAP average trace data
% load(fullfile(LOAD_PATH, filename_post), my_vars{:})
load(fullfile(path, filename), my_vars{:})
y_avg_sel = y_avg_cache_descending; % each cell entry is sorted stack of average traces for a frequency
y_avg_sel_cm = y_avg_cm_cache_descending;
%     y_avg_odd = cellfun(@plus, y_avg_sel, y_avg_sel_cm, 'UniformOutput', 0);
y_avg_odd = cellfun(@(x, y) (x + y)/2, y_avg_sel, y_avg_sel_cm, 'UniformOutput', 0);
%     y_avg_odd = cellfun(@minus, y_avg_sel, y_avg_sel_cm, 'UniformOutput', 0);
y_avg_even = cellfun(@(x, y) (x - y)/2, y_avg_sel, y_avg_sel_cm, 'UniformOutput', 0);

%     load(fullfile(path, filename_pre), 'y_avg_cache_descending')
%     y_avg_pre = y_avg_cache_descending;

% Stacked plot of cochlear microphonic
[fig, t] = tileplot_CAPABR(y_avg_sel_cm, all_freqs, A_descending);

title(t, [CHANNEL_KEY{i}, ' microphonic'], 'FontSize', 24)

t.TileSpacing = 'none';
t.Padding = 'tight';

% Maximize figure window size
fig.WindowState = 'maximized';    

% Save figure
% save_file = [filename_prestem, timepoint];
[~, save_file, ~] = fileparts(filename);
savefig_multiformat(gcf, SAVE_PATH, [save_file, '_', CHANNEL_KEY{i}, '_microphonic'])

%% Stacked plot of odd traces
[fig, t] = tileplot_CAPABR(y_avg_odd, all_freqs, A_descending);

title(t, [CHANNEL_KEY{i}, ' odd'], 'FontSize', 24)

t.TileSpacing = 'none';
t.Padding = 'tight';

% Maximize figure window size
fig.WindowState = 'maximized';    

% Save figure
% save_file = [filename_prestem, timepoint];
savefig_multiformat(gcf, SAVE_PATH, [save_file, '_', CHANNEL_KEY{i}, '_odd'])

%% Stacked plot of even traces
[fig, t] = tileplot_CAPABR(y_avg_even, all_freqs, A_descending);

title(t, [CHANNEL_KEY{i}, ' even'], 'FontSize', 24)

t.TileSpacing = 'none';
t.Padding = 'tight';

% Maximize figure window size
fig.WindowState = 'maximized';    

% Save figure
% save_file = [filename_prestem, timepoint];
savefig_multiformat(gcf, SAVE_PATH, [save_file, '_', CHANNEL_KEY{i}, '_even'])

end
    
    