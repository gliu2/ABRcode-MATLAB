% function plot_stack_humanABR(varargin)
% Load and plot a stack of averaged single trace human ABRs at different
% intensities, in ascending order of intensity.
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
FONT_SIZE = 18;
LINE_WIDTH = 1.5;
% DIRECTORY = 'D:\users\admin\Documents\George\Human Single Trial ABR\';
DIRECTORY = 'D:\George-abr\human_abr\';
PATH_X = 'd:\users\admin\Documents\GitHub\ABRcode-MATLAB\SINGS001_tiph_x_1-22-23.mat';
% ylim_max = [-0.5, 0.5];
% ylim_max = [-500, 500];
ylim_max = [-300, 300];
% ylim_max = [-120, 120];

disp('Select a folder with human ABR single trace files to load:')
path = uigetdir(DIRECTORY, 'Select folder with human ABR single trace files to analyze');
% [filename, path] = uigetfile(fullfile(DIRECTORY, '*.txt'));
load(PATH_X, 'x'); % import variable 'x' - single trace files do not contain metadata about x-axis values

% if nargin == 0
%     disp('Select a human ABR single trace file to load:')
%     [filename, path] = uigetfile(fullfile(DIRECTORY, '*.txt'));
%     load(PATH_X, 'x'); % import variable 'x' - single trace files do not contain metadata about x-axis values
% elseif nargin == 1
%     path_filename = varargin{1};
%     % Parse filename
%     [path, name, ext] = fileparts(path_filename);
%     filename = [name, ext];
%     load(PATH_X, 'x'); % import variable 'x' - single trace files do not contain metadata about x-axis values
% elseif nargin == 2
%     path_filename = varargin{1};
%     % Parse filename
%     [path, name, ext] = fileparts(path_filename);
%     filename = [name, ext];
%     x = varargin{2};
% end
% ind_underscore = strfind(filename, '_');
% filename_root = filename(1:ind_underscore(end));
% filename_tail = filename(end-6:end);
% filename_template = [filename_root, '\w*', filename_tail];

% Obtain directory listing of all ABR files in set (with same electrode
% setup, eg tip horizontal, tip vertical, or wick)
listing = dir(path);
fileList = {listing.name}';
num_files = length(fileList);
is_single_trace = zeros(num_files, 1);
is_rarefaction = zeros(num_files, 1);
is_condensation = zeros(num_files, 1);
is_vertical = zeros(num_files, 1);
db_levels = zeros(n_levels, 1);
% is_stack = zeros(num_files, 1);
for i=1:num_files
    filename = fileList{i};
    if length(filename) < 7
        continue
    end
    ind_underscore = strfind(filename, '_');
    filename_root = filename(1:ind_underscore(end));
    filename_root_last6 = filename_root(end-6:end);
    filename_tail = filename(end-6:end);
    
    is_single_trace(i) = strcmpi(filename_tail, '_st.TXT');
    is_rarefaction(i) = contains(filename_root_last6, 'R', IgnoreCase=true);
    is_condensation(i) = contains(filename_root_last6, 'C', IgnoreCase=true);
    is_vertical(i) = strcmpi(filename_root(end-1), 'v');
    
%     is_stack(i) = any(regexp(fileList{i}, filename_template));
end

% Obtain decibel level of files
n_levels = sum(is_stack);
db_levels = zeros(n_levels, 1);
ind_stack = find(is_stack);
stack = cell(n_levels, 1);
for i=1:n_levels
    filename = fileList{ind_stack(i)};
    nameNumber = filename(ind_underscore(end) + 1: end - 7);
    db_levels(i, 1) = nameToNumber(nameNumber);
    stack{i} = filename;
end

% Sort files in ascending order
[db_levels_descend, ind_descend] = sort(db_levels, 'descend'); % default is ascending order
stack_descend = stack(ind_descend);

%% Plot stack of averaged single trace ABRs in ascending order
dB = cell(n_levels, 1);
dB(:) = {'dB'};
A_descending_cell = cellstr(num2str(db_levels_descend));
newYlabels=cellfun(@(x,y) [x  ' ' y], A_descending_cell, dB, 'un', 0);
this_xlim = [];

fig = figure;
t = tiledlayout(n_levels, 1, 'TileIndexing', 'columnmajor');

for i=1:n_levels
    filename = fullfile(path, stack_descend{i});
%     A = load_humanABR_singletracedata(this_filename);
    A = load_humanABR_singletracedata(filename, HIGHPASS_FREQ, LOWPASS_FREQ, NOTCHFILTER_FREQ);
    y = mean(A, 2);
    y = y*CONVERSION_SINGLETRACE_NV; % convert to units of uV
    
    nexttile(t)
    plot(x, y, 'LineWidth', 3)
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
    
    % Mark wave 1 peak and following trough
    if i==1
        % Get wave 1 peak and following trough at highest stimulus
        % level
        [peak_pt, trough_pt, amp, lat_peak, lat_trough] = get_wave1_humanaverageABR(x, y);
    else
        % Get wave 1 peak and following trough at lower stimulus
        % level. Ensure peak latency does not decrease.
        [peak_pt, trough_pt, amp, lat_peak, lat_trough] = get_wave1_humanaverageABR(x, y, peak_pt(1), trough_pt(1));
    end

    hold on
    CIRCLE_SIZE = 54; % default 36
    scatter(peak_pt(1), peak_pt(2), CIRCLE_SIZE, 'green', 'filled') % peak
    scatter(trough_pt(1), trough_pt(2), CIRCLE_SIZE, 'red', 'filled') % trough
    hold off

    % Display wave 1 measurements in corner of plot
    xL=xlim;
    yL=ylim;
    str = sprintf('P-P_I = %.0f nV', amp);
    text(0.99*xL(2), 0.99*yL(1), str,'HorizontalAlignment','right','VerticalAlignment','bottom')

    ylabel(newYlabels{i}, 'FontSize', FONT_SIZE)
    % rotate y label to make horizontal
    hYLabel = get(gca,'YLabel');
    set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right', 'FontSize', FONT_SIZE)

    % Remove extraneous axis tick marks and x-axis from all but bottom
    % tile
    set(gca,'box','off')
    if i~=n_levels
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'XColor','none')
    end

    % Set text size of current axis
    set(gca,'FontSize',24)
end
     
xlabel('Time (ms)', 'FontSize', FONT_SIZE)
 
t.TileSpacing = 'none';
t.Padding = 'tight';
plot_title = filename_root(1:end-1);
title(t, plot_title, 'FontSize', FONT_SIZE, 'Interpreter','none')
ylabel(t, 'Amplitude (nV)', 'FontSize', FONT_SIZE)


% % Maximize figure window size
% fig.WindowState = 'maximized';

% % Save figure
% %     disp('Saving figure')
% [~, save_file, ~] = fileparts(filename);
% savefig_multiformat(gcf, SAVE_PATH, save_file)


% end