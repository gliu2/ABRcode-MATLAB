% function plot_stack_humanABR_fullstack(varargin)
% Load and plot a stack of averaged single trace human ABRs at different
% intensities, in ascending order of intensity.
%
% Only analyzes S001 data. Cell 2 uses syntax of filenames to input single
% trace human ABR data for S001, and performs offline filtering and
% plotting to create figure. See other .m file for analysis of other
% subject data which has different file formatting. 
%
% Optional input 1: path_filename (full path and filename to block average
% / single trace file.
%
% Optional input 2: X-axis time values (ms) 
%
% Last edit George Liu 8-31-23
%
% Dependencies: load_humanABR_average.m, load_humanABR_singletracedata.m,
% nameToNumber.m

% CONVERSION_SINGLETRACE_UV = 0.0030; % multiply single trace data by this factor to put in units of uV
CONVERSION_SINGLETRACE_NV = 3; % multiply single trace data by this factor to put in units of nV
HIGHPASS_FREQ = 300; % Hz
LOWPASS_FREQ = 3000; % Hz
NOTCHFILTER_FREQ = 60; % Hz
FONT_SIZE = 24;
LINE_WIDTH = 1.5;
% DIRECTORY = 'D:\users\admin\Documents\George\Human Single Trial ABR\';
SAVE_PATH = 'd:\users\admin\Documents\George\Results_analyzed'; % path for saving figures
DIRECTORY = 'D:\George-abr\human_abr\';
PATH_X = 'd:\users\admin\Documents\GitHub\ABRcode-MATLAB\SINGS001_tiph_x_1-22-23.mat';
% ylim_max = [-0.5, 0.5];
% ylim_max = [-500, 500];
ylim_max = [-300, 300];
% ylim_max = [-120, 120];
myColors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980]}; % blue and orange


% SINGLE TRACE KEYWORD
SINGLETRACE_KEY = 'blk'; % bulk average, i.e. average of rarefaction and condensation traces

%% Hard code these parameters
% MONTAGE = 'tiph';
MONTAGE = 'tipv';
% MONTAGE = 'wick';
path = 'D:\George-abr\human_abr\SingleTrials_aas_01\Matt Initial Recordings';
load(PATH_X, 'x'); % import variable 'x' - single trace files do not contain metadata about x-axis values

% Obtain directory listing of all ABR files in set (with same electrode
% setup, eg tip horizontal, tip vertical, or wick)
listing = dir(path);
fileList = {listing.name}';
num_files = length(fileList);
is_single_trace = zeros(num_files, 1);
is_montage = zeros(num_files, 1);
is_stack = zeros(num_files, 1);
for i=1:num_files
    filename = fileList{i};
    is_single_trace(i) = contains(filename, SINGLETRACE_KEY, IgnoreCase=true);
    is_montage(i) = contains(filename, MONTAGE, IgnoreCase=true);
    is_stack(i) = is_single_trace(i) && is_montage(i);
end

% Obtain decibel level of files
n_levels = sum(is_stack);
db_levels = zeros(n_levels, 1);
ind_stack = find(is_stack);
stack = cell(n_levels, 1);
for i=1:n_levels
    filename = fileList{ind_stack(i)};
    ind_underscore = strfind(filename, '_');
    nameNumber = filename(ind_underscore(end) + 1: end - 7);
    db_levels(i, 1) = nameToNumber(nameNumber);
    stack{i} = filename;
end

% Sort files in descending order
[db_levels_descend, ind_descend] = sort(db_levels, 'descend'); % default is ascending order
stack_descend = stack(ind_descend);

%% Plot stack of averaged single trace ABRs in descending order
dB = cell(n_levels, 1);
dB(:) = {'dB (nV)'};
A_descending_cell = cellstr(num2str(db_levels_descend));
newYlabels=cellfun(@(x,y) [x  ' ' y], A_descending_cell, dB, 'un', 0);
this_xlim = [];
X_csv_cache = cell(n_levels, 1);
y_avg_cache = cell(n_levels, 1);

fig = figure;
% t = tiledlayout(n_levels, 1, 'TileIndexing', 'columnmajor');
t = tiledlayout(n_levels - 1, 1, 'TileIndexing', 'columnmajor'); % exclude -30 dB

% for i=1:n_levels
for i=1:n_levels - 1 % exclude -30 dB
    filename = fullfile(path, stack_descend{i});
    A = load_humanABR_singletracedata(filename);
%     A = load_humanABR_singletracedata(filename, HIGHPASS_FREQ, LOWPASS_FREQ, NOTCHFILTER_FREQ);
    
    % Filter trace
    N_ORDER = 700; % Half mainlobe width (-3 dB, approximation of cutoff frequency) of 226.484 Hz
%     A = filter_highpass_blackman(A, N_ORDER); % Blackman high pass filter
%     A = filter_highpass_blackman(A, 520); % Blackman high pass filter cutoff 300 Hz
%     A = filter_highpass_blackman(A, 300); % Blackman high pass filter cutoff 500 Hz
    A = filter_highpass_blackman(A, 200); % Blackman high pass filter cutoff 750 Hz
    A = filter_lowpass_blackman(A, 20); % length 20 gives Low pass frequency cutoff of 8.4 kHz
%     A = filter_lowpass_blackman(A, 26); % length 26 gives Low pass frequency cutoff of 6.4 kHz
%     A = filter_lowpass_blackman(A, 32); % length 32 gives Low pass frequency cutoff of 5.1 kHz
%     A = filter_lowpass_blackman(A, 40); % length 40 gives Low pass frequency cutoff of 4.1 kHz
%     A = filter_lowpass_blackman(A, 54); % length 54 gives Low pass frequency cutoff of 3 kHz
%     A = filter_lowpass_blackman(A, 150); % Blackman low pass filter with cutoff of 1.05 kHz (Blackman window of order 150)
%     A = filter_lowpass_blackman(A, 150);  % Apply twice to suppress cycle (and alias) noise more

    % Keep only portion of A between 0 and 10 ms
    x_range = [0, 10]; % ms time window to analyze
    is_x_inrange = x >= x_range(1) & x <= x_range(2);
    
    x_inrange = x(is_x_inrange);
    A = A(is_x_inrange, :);

    y = mean(A, 2);
    y = y*CONVERSION_SINGLETRACE_NV; % convert to units of uV
    
    nexttile(t)
%     plot(x, y, 'LineWidth', 3) % Plot full time window including before stimulus onset
    plot(x_inrange, y, 'LineWidth', 3)
    hold on
    xline(0, '--', 'LineWidth', LINE_WIDTH)
    hold off
    ylim(ylim_max)
    
%     % Set x axis limits to match limits of highest stimulus level (1st row)
%     if i==1
%         this_xlim = xlim;
%     else
%         xlim(this_xlim);
%     end
    
%     % Mark wave 1 peak and following trough
%     if i==1
%         % Get wave 1 peak and following trough at highest stimulus
%         % level
%         [peak_pt, trough_pt, amp, lat_peak, lat_trough] = get_wave1_humanaverageABR(x, y);
%     else
%         % Get wave 1 peak and following trough at lower stimulus
%         % level. Ensure peak latency does not decrease.
%         [peak_pt, trough_pt, amp, lat_peak, lat_trough] = get_wave1_humanaverageABR(x, y, peak_pt(1), trough_pt(1));
%     end
% 
%     hold on
%     CIRCLE_SIZE = 54; % default 36
%     scatter(peak_pt(1), peak_pt(2), CIRCLE_SIZE, 'green', 'filled') % peak
%     scatter(trough_pt(1), trough_pt(2), CIRCLE_SIZE, 'red', 'filled') % trough
%     hold off
% 
%     % Display wave 1 measurements in corner of plot
%     xL=xlim;
%     yL=ylim;
%     str = sprintf('P-P_I = %.0f nV', amp);
%     text(0.99*xL(2), 0.99*yL(1), str,'HorizontalAlignment','right','VerticalAlignment','bottom')
% 
    ylabel(newYlabels{i}, 'FontSize', FONT_SIZE)
    % rotate y label to make horizontal
    hYLabel = get(gca,'YLabel');
    set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right', 'FontSize', FONT_SIZE)

    % Remove extraneous axis tick marks and x-axis from all but bottom
    % tile
    set(gca,'box','off')
%     if i~=n_levels
    if i~=n_levels - 1
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'XColor','none')
    end

    % Set text size of current axis
    set(gca,'FontSize',24)
    
    % Cache variables
    X_csv_cache{i} = A;
    y_avg_cache{i} = y;
end
     
xlabel('Time (ms)', 'FontSize', FONT_SIZE)
 
t.TileSpacing = 'none';
t.Padding = 'tight';
[~, filename_noext, ~] = fileparts(filename);
ind_underscore_this = strfind(filename_noext, '_');
plot_title = filename_noext(1:ind_underscore_this (end)-1);
title(t, plot_title, 'FontSize', FONT_SIZE, 'Interpreter','none')
ylabel(t, 'Stimulus level (dB HL)', 'FontSize', FONT_SIZE)


% Maximize figure window size
% fig.WindowState = 'maximized';

% Save figure
%     disp('Saving figure')
% [~, save_file, ~] = fileparts(filename);
save_file = [plot_title, '_humanABRstackfull'];
savefig_multiformat(gcf, SAVE_PATH, save_file)

%% Analyze threshold
cutoff = 1;
[threshold_D, metric_D, dist_innerprod] = get_thresh_CAPABR_D(X_csv_cache(1:end-1), db_levels_descend(1:end-1), cutoff, true);
xlabel('Stimulus level (dB HL')

%% Save threshold plot
save_file2 = [plot_title, '_humanABRstackfull_threshold'];
savefig_multiformat(gcf, SAVE_PATH, save_file2)

%% Fig 3: Plot human average traces in tiled plot
fig2 = figure;
NUM_LEVELS = n_levels - 1;
A_descending = db_levels_descend(1:end-1);
N_COLUMNS = 3;
t2 = tiledlayout(NUM_LEVELS, N_COLUMNS, 'TileIndexing', 'columnmajor');

% Make sure bounds of plot are within y limits
ymax_data = max(cell2mat(y_avg_cache), [], "all");
ymin_data = min(cell2mat(y_avg_cache), [], "all");
ylim_max = [ymin_data, ymax_data];

for j = 1:NUM_LEVELS
    y_avg = y_avg_cache{j};
    A_csv = db_levels_descend(j);
    
    nexttile(t2)
%     plot(X, y_avg, 'LineWidth', 3)
    stdshade(X_csv_cache{j}', 0.1, myColors{1}, x_inrange, [], 3);
%     steshade(X_csv_cache{j}', 0.1, myColors{1}, X, [], 3);
    ylim(ylim_max)

    % show ylabels for first column only
    ylabel([num2str(A_csv), ' dB (nV)'], 'FontSize', 24)
    % rotate y label to make horizontal
    hYLabel = get(gca,'YLabel');
    set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right', 'FontSize', 24)
    
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
    set(gca,'FontSize',24)

    % Show title of column one
    if mod(j, NUM_LEVELS)==1
        title(['Human ABR']) 
    end
    
    xlabel('Time (ms)', 'FontSize', 24)
    xlim([0, 13])
%     xlim([0, 10])
    ylim([-375, 400])
end

% Plot distributions of inner products
max_innerprod = max(dist_innerprod, [], "all");
min_innerprod = min(dist_innerprod, [], "all");
innerprod_xlim = [min_innerprod, max_innerprod];
max_bin_val = zeros(NUM_LEVELS, 1);
axesHandle = zeros(NUM_LEVELS, 1);

LINE_WIDTH = 2;
for j = 1:NUM_LEVELS
    axesHandle(j) = nexttile(t2);
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

    xlabel('Inner product (nV^2)')

    % plot mean as blue line
%     line([signal_mean signal_mean], ylim, 'LineWidth', 1, 'Color', 'b'); % vertical line for cutoff
    hold off

    if j==1
        title('Histogram') 
    end

    % show ylabels for first column only
    if j==5
        ylabel('Density', 'FontSize', 24)
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
    set(gca,'FontSize',24)
    
    % Same x axes for all histograms
%     xlim(innerprod_xlim)
    xlim([-2, 3]*10^6)
    max_bin_val(j) = max(h.Values, [], "all");
    
    % Set y axis to log
    set(gca,'YScale','log')
end
% Set y axes same
set(axesHandle, 'YLim', [0, max(max_bin_val)]);

% Plot D'
nexttile([NUM_LEVELS, 1])
plot(metric_D, A_descending, 'LineWidth', 3)
xlabel('Cohen D (a.u.)', 'FontSize', 24)
my_xline = xline(cutoff, '--', ['Cutoff @ D = ', num2str(cutoff)], 'LineWidth', LINE_WIDTH);
% title(['Threshold at ', num2str(round(threshold_D, 1)), ' dB'], 'FontSize', 24)
title(['Threshold'], 'FontSize', 24)
ylabel('Stimulus level (dB HL)', 'FontSize', 24)
set(my_xline, 'FontSize', 18)
set(gca,'TickDir','out');
set(gca,'box','off')
% Set text size of current axis
set(gca,'FontSize',24)
ylim([A_descending(end)-5, A_descending(1)+5]);
xlim_bounds = [min([metric_D; 0]), 1.05*max(metric_D)];
xlim(xlim_bounds);

% Draw horizontal line at threshold amplitude 
if ~isnan(threshold_D)
    hold on
    my_yline = yline(threshold_D, '--', ['Threshold = ', num2str(round(threshold_D, 1)), ' dB HL'], 'LineWidth', LINE_WIDTH, 'Color', 'k'); % horizontal line at exact threshold
    set(my_yline, 'FontSize', 18)
    hold off
else
    disp('Warning: No ABR threshold detected!')
end

% Format figure
% xlabel(t, 'Time (ms)', 'FontSize', 24)
ylabel(t2, 'Stimulus level (dB HL)', 'FontSize', 24)
% title(t, CHANNEL_KEY{i}, 'FontSize', 24)

t2.TileSpacing = 'none';
t2.Padding = 'tight';

% Maximize figure window size
fig2.WindowState = 'maximized';    

%% Save figure
save_file = [plot_title, '_humanABRstackfull'];
savefig_multiformat(gcf, SAVE_PATH, [save_file, '_thresholdpaper_Fig3'])


% end