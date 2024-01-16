% function plot_threshold_shift(varargin)
% Plot threshold shift for one mouse
%
% 7/6/2023 George Liu
% Dependencies: savefig_multiformat.m

SAVE_FIGURES = 1; % Hard code logical for whether or not to save figures 

MAX_TIMEPOINTS = 4;
% LEG_LABEL = ["Pre", "10 min", "30 min", "60 min"];
CHANNEL_KEY = {'ABR', 'CAP'};
LOAD_PATH = 'd:\users\admin\Documents\George\Results';
SAVE_PATH = 'd:\users\admin\Documents\George\Results_analyzed'; % path for saving figures

%%
[filename, path] = uigetfile(fullfile(LOAD_PATH, '*.mat')); % returns 0 for outputs if no selection
display(['Opening ', fullfile(path, filename)])
% Analyze metadata in single trial ABR/CAP filename
channel = nan;
for i=1:length(CHANNEL_KEY)
    if any(strfind(filename, CHANNEL_KEY{i}))
        channel = i;
    end
end

underscore_loc = strfind(filename, '_');

filename_prestem = filename(1 : underscore_loc(2) - 1);
listing = dir(fullfile(path, [filename_prestem, '*', CHANNEL_KEY{channel}, '_results.mat']));

allnames = {listing.name}'; % column cell array of filenames

%%
n_timepoints = length(allnames);

n_thresholds = cell(n_timepoints, 1);
n_metrics = cell(n_timepoints, 1);
timepoint = nan(n_timepoints, 1);
leg_timepoint = strings(n_timepoints,1);
for i = 1:n_timepoints
    this_filename = allnames{i};
    this_fullfilename = fullfile(path, this_filename);
    display(['Opening ', this_fullfilename])
    load(this_fullfilename, "A_descending", "all_freqs", "metric", "threshold")
    
    n_thresholds{i} = threshold;
    n_metrics{i} = metric;
    
    % Get metadata to determine time point
    underscore_loc = strfind(this_filename, '_');
    filename_time = this_filename(underscore_loc(2) + 1 : underscore_loc(3) - 1);
    dash_loc = strfind(filename_time, '-');
    if any(dash_loc)
        filename_time = filename_time(1 : dash_loc(1) - 1);
    end
    
    if any(strfind(filename_time, "pre"))
        timepoint(i) = 0;
        leg_timepoint(i) = "Pre";
    elseif any(strfind(filename_time, "10"))
        timepoint(i) = 10;
        leg_timepoint(i) = "10 min";
    elseif any(strfind(filename_time, "30"))
        timepoint(i) = 30;
        leg_timepoint(i) = "30 min";
    elseif any(strfind(filename_time, "60"))
        timepoint(i) = 60;
        leg_timepoint(i) = "60 min";
    end
end

% Sort plots in order from earlier to later
[~, I] = sort(timepoint);
metrics_sorted = n_metrics(I);
thresholds_sorted = n_thresholds(I);
leg_timepoint_sorted = leg_timepoint(I);
% metrics_sorted = n_metrics(I);

%%
XLIM = [7, 35];
YLIM = [0, 100];

all_freqs_khz = all_freqs/1000;

% % Get new colormap with same hue, different saturation for post-injection
% % thresholds
% originalrgb = [1 0 1]; %replace by whatever rgb colour you want
% originalhsv = rgb2hsv(originalrgb);  %get the HSV values of your original colour. We really only care about the hue
% maphsv = rgb2hsv(gray);  %without any argument gray returns 64 values; convert to hsv
% maphsv(:, 1) = originalhsv(1);  %replace gray hue by original hue
% maphsv(:, 2) = originalhsv(2); %replace saturation. Anything but 0 will work
% newmap = hsv2rgb(maphsv);
% % colormap(newmap);
% 
% 
% colormap_linspace = floor(linspace(1, size(newmap, 1), n_timepoints));
% colors = newmap(colormap_linspace, :);
% % colors = [[0, 0, 0]; newmap(colormap_linspace(2:end), :)];

colors = [[0.5, 0.5, 0.5]; flipud(copper(n_timepoints-1))];

fig = figure;
hold on
for i=1:n_timepoints
    semilogx(all_freqs_khz, thresholds_sorted{i}, 'LineWidth', 3, 'Color', colors(i, :))
end
hold off
legend(leg_timepoint_sorted, 'Location', 'southeast') 

set(gca, 'XScale', 'log');
xticks(all_freqs_khz)
xlim(XLIM)
ylim(YLIM)


xlabel('Frequency (kHz)', 'FontSize', 24)
ylabel('Threshold (dB)', 'FontSize', 24)
% title(t, CHANNEL_KEY{i}, 'FontSize', 24)

title([filename_prestem,  '_', CHANNEL_KEY{channel}], 'FontSize', 24, 'Interpreter', 'none')

% Set text size of current axis
set(gca,'FontSize',24)
            

% Maximize figure window size
fig.WindowState = 'maximized';

if SAVE_FIGURES
    % Save figure
    %     disp('Saving figure')
    [~, save_file, ~] = fileparts(filename_prestem);
    savefig_multiformat(gcf, SAVE_PATH, [save_file, '_', CHANNEL_KEY{channel}, '_thresholdshift'])
end

%% Save data
savefilename = [filename_prestem, '_thresholdshift.mat'];
save(fullfile(SAVE_PATH, savefilename), 'all_freqs_khz', 'thresholds_sorted', 'leg_timepoint_sorted', 'metrics_sorted')