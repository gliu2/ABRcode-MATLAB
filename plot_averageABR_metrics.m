% plot_averageABR_metrics
% script 
%
% **IMPORTANT** BEFORE RUNNING MANUALLY SET THE FOLLOWING VARIABLES:
% IS_NOISE
% ANALYZE_FULLDATA_MICE_ONLY
% IS_RIGHTEARONLY
%
% 8/2/21 George Liu

%% MANUALLY SET THESE SWITCHES
IS_NOISE = 1;
ANALYZE_FULLDATA_MICE_ONLY = 0;
IS_RIGHTEARONLY = 0;
IS_LEFTEARONLY = 0;

%% Constants
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures
VARNAMES = {'Date', 'Amplitude_8khz', 'Latency_8khz', 'ABRThreshold_8khz', 'DPOAEThreshold_8khz', 'Amplitude_16khz', 'Latency_16khz', 'ABRThreshold_16khz', 'DPOAEThreshold_16khz', 'Amplitude_32khz', 'Latency_32khz', 'ABRThreshold_32khz', 'DPOAEThreshold_32khz'};
FREQUENCIES = [8, 16, 32];
NAME_ABRTHRESH = {'ABRThreshold_8khz', 'ABRThreshold_16khz', 'ABRThreshold_32khz'};
NAME_DPOAETHRESH = {'DPOAEThreshold_8khz', 'DPOAEThreshold_16khz', 'DPOAEThreshold_32khz'};
NAME_AMPLITUDE = {'Amplitude_8khz', 'Amplitude_16khz', 'Amplitude_32khz'};
NAME_LATENCY = {'Latency_8khz', 'Latency_16khz', 'Latency_32khz'};
METRICS = {NAME_ABRTHRESH, NAME_DPOAETHRESH, NAME_AMPLITUDE, NAME_LATENCY};
METRIC_TITLES = {'ABR Threshold', 'DPOAE Threshold', 'Wave 1 amplitude', 'Wave 1 latency'};
MEAN_SE_STRING = ' (mean +/- SE)';
METRIC_YLABELS = {'Threshold (dB)', 'Threshold (dB)', 'Amplitude (nV)', 'Latency (ms)'};
MY_XLIM = [7, 36];
METRIC_YLIMS = {[20, 90], [20, 90], [0, 5100], [1.3, 1.7]};
METRIC_LEGENDLOC = {'northwest', 'southwest', 'southwest', 'southwest'};
NAME_TIMES = {'Pre', '1 day', '1 week'};
X_TIMES = [-3, 1, 7];

%% Select excel file with aggregate average ABR analysis results
[filename,path] = uigetfile('d:\users\admin\Documents\George\aggregate_results\*.xlsx',...
    'Select an Excel File');

% Load data
T = readtable(fullfile(path, filename),...
    'Range','A3',...
    'ReadVariableNames',true);

T_key = T(:, 1:4);
T_pre = T(:, 5:17);
T_24h = T(:, 19:31);
T_1w = T(:, 33:end);

T_pre.Properties.VariableNames = VARNAMES;
T_24h.Properties.VariableNames = VARNAMES;
T_1w.Properties.VariableNames = VARNAMES;

%% Optional: analyze only mice with full set of data at all 3 time points
if ANALYZE_FULLDATA_MICE_ONLY
    ismouse_fulldata = ~any(isnan(T_pre{:,:}), 2) & ~any(isnan(T_24h{:,:}), 2) & ~any(isnan(T_1w{:,:}), 2);
else
    ismouse_fulldata = ones(size(T_pre{:,:}, 1), 1);
end

if IS_RIGHTEARONLY
    ismouse_fulldata = ismouse_fulldata & contains(T_key.Side, 'right');
elseif IS_LEFTEARONLY
    ismouse_fulldata = ismouse_fulldata & contains(T_key.Side, 'left');
end
ind_ismouse_fulldata = find(ismouse_fulldata); % Convert logical indexing to row coords

%% Plot metrics vs frequency, one plot per time point

num_metrics = length(METRICS);
figure
t = tiledlayout(2, 2);

% Cycle through each metric
for j=1:num_metrics
    this_metric = METRICS{j};
    
    % For this metric, plot 3 curves (pre, 1 day, and 1 week)
    T_all = {T_pre, T_24h, T_1w};
    num_timepoints = length(T_all);
    counts_cache = zeros(num_timepoints,3);
    for i=1:num_timepoints
        this_T = T_all{i};
        T_metric = this_T(ind_ismouse_fulldata, this_metric);
        metric_mean = mean(T_metric{:,:}, 1, 'omitnan');
        metric_std = std(T_metric{:,:}, 1, 'omitnan');
        counts = sum(~isnan(T_metric{:,:}), 1);
        metric_ste = metric_std./sqrt(counts);

        if i==1
            nexttile
        elseif i==2
            hold on
        end

        errorbar(FREQUENCIES, metric_mean, metric_ste)

        if i==num_timepoints
            hold off
        end

        counts_cache(i,:) = counts;

    end

    hAx = gca;
    hAx.XScale='log';
    xlim(MY_XLIM)
    ylim(METRIC_YLIMS{j});
    xticks(FREQUENCIES)
    xticklabels({FREQUENCIES})
%     xlabel('Frequency (kHz)')
    ylabel(METRIC_YLABELS{j})
    title([METRIC_TITLES{j}, MEAN_SE_STRING])

    counts_mean = num2cell(round(mean(counts_cache, 2))');
    newLegendlabels=cellfun(@(x,y) [x,  ' (n=', num2str(y), ')'], NAME_TIMES, counts_mean, 'un', 0);
    legend(newLegendlabels, 'Location', METRIC_LEGENDLOC{j})
end

t.TileSpacing = 'tight';
t.Padding = 'tight';
xlabel(t, 'Frequency (kHz)')
if IS_NOISE
    title(t, 'Noise exposed')
else
    title(t, 'Control')
end

%% Plot individual mice data
% Amplitudes
this_metric = METRICS{3};
metric_pre = T_pre(ind_ismouse_fulldata, this_metric);
metric_24h = T_24h(ind_ismouse_fulldata, this_metric);
metric_1w = T_1w(ind_ismouse_fulldata, this_metric);

delta_24h = metric_24h{:,:} - metric_pre{:,:};
delta_1w = metric_1w{:,:} - metric_pre{:,:};

counts = sum(~isnan(metric_pre{ind_ismouse_fulldata,:}), 1);
these_mice = T_key.Mouse(ind_ismouse_fulldata);
num_these_mice = length(unique(these_mice));
% % Add random jitter
% rand_jitter = 300*rand(length(T.Frequency),1);
% group_slide = [zeros(9,1); ones(9,1)]*1000;

figure
t2 = tiledlayout(1, 3);
for ff=1:3
    nexttile
    
    cache_names = cell(counts(ff), 1);
    
    % create a color map for mice
    num_colors = num_these_mice;
    colors_p = jet(num_colors);
    for k=1:counts(ff)
        this_mouse = these_mice{k};
        ind_this_mouse = find(contains(unique(these_mice), this_mouse));
        y = [metric_pre{k,ff}, metric_24h{k,ff}, metric_1w{k,ff}];
        if contains(T_key.Side(ind_ismouse_fulldata(k)), 'right')
            plot(X_TIMES, y, '-o', 'color', colors_p(ind_this_mouse,:))
        else
            plot(X_TIMES, y, '-x', 'color', colors_p(ind_this_mouse,:))
        end
        if k==1
            hold on
        end
        
        cache_names{k} = this_mouse;
    end
    hold off
    ylim([0, 8200])
    title([num2str(FREQUENCIES(ff)), ' kHz'])
    xticks(X_TIMES)
    xticklabels(NAME_TIMES)
    legend(cache_names, 'Location', 'northeast')
end
xlabel(t2, 'Time after exposure')
ylabel(t2, METRIC_YLABELS{3})
if IS_NOISE
    this_title = [METRIC_TITLES{3}, ', noise exposed'];
else
    this_title = [METRIC_TITLES{3}, ', control'];
end
title(t2, this_title)

