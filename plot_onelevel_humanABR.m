function plot_onelevel_humanABR(varargin)
% Load and plot a stack of averaged single trace human ABRs at different
% intensities, in ascending order of intensity.
%
% Optional input 1: path_filename (full path and filename to block average
% / single trace file.
%
% Optional input 2: X-axis time values (ms) 
%
% Last edit George Liu 2-26-23
%
% Dependencies: load_humanABR_singletracedata.m, steshade.m 

% CONVERSION_SINGLETRACE_UV = 0.0030; % multiply single trace data by this factor to put in units of uV
CONVERSION_SINGLETRACE_NV = 3; % multiply single trace data by this factor to put in units of nV
HIGHPASS_FREQ = 300; % Hz
LOWPASS_FREQ = 3000; % Hz
NOTCHFILTER_FREQ = 60; % Hz
FONT_SIZE = 16;
LINE_WIDTH = 1.5;
% DIRECTORY = 'D:\users\admin\Documents\George\Human Single Trial ABR\';
DIRECTORY = 'D:\George-abr\human_abr\';
PATH_X = 'd:\users\admin\Documents\GitHub\ABRcode-MATLAB\SINGS001_tiph_x_1-22-23.mat';
% ylim_max = [-0.5, 0.5];
% ylim_max = [-500, 500];
ylim_max = [-1500, 1500];
xlim_max = [-2, 12];
% ylim_max = [-120, 120];

% Select rarefaction or condensation single trace ABR file
if nargin == 0
    disp('Select a human ABR single trace file to load:')
    [filename, path] = uigetfile(fullfile(DIRECTORY, '*.txt'));
    load(PATH_X, 'x'); % import variable 'x' - single trace files do not contain metadata about x-axis values
elseif nargin == 1
    path_filename = varargin{1};
    % Parse filename
    [path, name, ext] = fileparts(path_filename);
    filename = [name, ext];
    load(PATH_X, 'x'); % import variable 'x' - single trace files do not contain metadata about x-axis values
elseif nargin == 2
    path_filename = varargin{1};
    % Parse filename
    [path, name, ext] = fileparts(path_filename);
    filename = [name, ext];
    x = varargin{2};
end

% Select condensation or rarefaction single trace ABR file (whichever was
% not selected before)
disp('Select a human ABR single trace file to load (condensation if prev selected was rarefaction, or vice versa):')
[filename2, path2] = uigetfile(fullfile(DIRECTORY, '*.txt'));

% Average rarefaction and condensation traces
A1 = load_humanABR_singletracedata(fullfile(path, filename), HIGHPASS_FREQ, LOWPASS_FREQ, NOTCHFILTER_FREQ);
A2 = load_humanABR_singletracedata(fullfile(path2, filename2), HIGHPASS_FREQ, LOWPASS_FREQ, NOTCHFILTER_FREQ);
A_avg = (A1 + A2)/2;

% Plot average
A_avg = A_avg*CONVERSION_SINGLETRACE_NV; % convert to units of nV
figure, stdshade(A_avg', 0.3, 'b', x, []);


xlim(xlim_max)
ylim(ylim_max)
set(gca,'FontSize', FONT_SIZE)
ylabel('Amplitude (nV)')
xlabel('Time (ms)')

% Mark wave 1 peak and following trough
y = mean(A_avg, 2);
[peak_pt, trough_pt, amp, lat_peak, lat_trough] = get_wave1_humanaverageABR(x, y);

hold on
CIRCLE_SIZE = 54; % default 36
scatter(peak_pt(1), peak_pt(2), CIRCLE_SIZE, 'green', 'filled') % peak
scatter(trough_pt(1), trough_pt(2), CIRCLE_SIZE, 'red', 'filled') % trough
hold off

% Obtain wave 1 standard deviation for single traces' wave 1 amps
% y_std = std(A_avg, [], 2);
% wave1std = y_std(x==lat_peak);
wave1peak_dist = A_avg(x==lat_peak, :);
wave1trough_dist = A_avg(x==lat_trough, :);
wave1amp_dist = wave1peak_dist - wave1trough_dist;
wave1std = std(wave1amp_dist);

% Display wave 1 measurements in corner of plot
xL=xlim;
yL=ylim;
str = sprintf('Wave I amp = %.0f nV', amp);
str2 = sprintf('Wave I lat = %.2f ms', lat_peak);
str3 = sprintf('Wave I std = %.0f nV', wave1std);
text(0.99*xL(2), 0.79*yL(1), str,'HorizontalAlignment','right','VerticalAlignment','bottom')
text(0.99*xL(2), 0.89*yL(1), str2,'HorizontalAlignment','right','VerticalAlignment','bottom')
text(0.99*xL(2), 0.99*yL(1), str3,'HorizontalAlignment','right','VerticalAlignment','bottom')

% Legend
legend({'STD', 'mean', 'I peak', 'I trough'})
legend('boxoff')

end