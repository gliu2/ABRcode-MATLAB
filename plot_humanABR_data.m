% plot_humanABR_data(varargin)
%Script for loading and visualizing average trace human ABR data
%
% Last edit George Liu 2-11-23
% Dependencies: load_humanABR_average.m, load_humanABR_singletracedata.m

FILE_NAME_SINGLE = 'D:\users\admin\Documents\George\Human Single Trial ABR\SINGS001_tiph_80blk.TXT';
FILE_NAME_AVERAGE = 'D:\users\admin\Documents\George\Human Single Trial ABR\SINGS001_tiph_80avg.TXT';
TIME_UNIT = 'ms';
VOLTAGE_UNIT = 'uV';
CONVERSION_SINGLETRACE_UV = 0.0030; % multiply single trace data by this factor to put in units of uV
FONT_SIZE = 18;
LINE_WIDTH = 1.5;

% Load average data
[x, y, y_uv] = load_humanABR_average(FILE_NAME_AVERAGE);

% Load single trace data
A = load_humanABR_singletracedata(FILE_NAME_SINGLE); % 1024 x 500 matrix, 500 single traces each with 1024 samples
A = A * CONVERSION_SINGLETRACE_UV;

A_avg = mean(A, 2); 

% uv_conversion = y_uv(1) / A_avg(1);
% A_avg_uV = A_avg * CONVERSION_SINGLETRACE_UV;
A_avg_uV = A_avg;

% Plot average data
f1 = figure;
set(gca, 'FontSize', FONT_SIZE)
plot(x, y_uv, 'b-', 'LineWidth', LINE_WIDTH)
hold on
plot(x, A_avg_uV, 'r-', 'LineWidth', LINE_WIDTH)
xline(0, '--', 'LineWidth', LINE_WIDTH)
hold off

xlabel('ms', 'FontSize', FONT_SIZE)
ylabel('Amplitude (uV)', 'FontSize', FONT_SIZE)
legend({'Average trace', 'Single traces'}, 'FontSize', FONT_SIZE, 'location', 'SouthWest')
title('S001, tip-h, 80 dB, n=500 single traces', 'FontSize', FONT_SIZE)

% Plot single trace data
f2 = figure;
n_traces = size(A, 2);
cmap = summer(n_traces);
hold on
for k=1:n_traces
    plot(x, A(:, k), 'Color', [cmap(k, 1), cmap(k, 2), cmap(k, 3), 0.4]) % transparency set to 4th color value
end
plot(x, A_avg_uV, 'r-', 'LineWidth', LINE_WIDTH)
xline(0, '--', 'LineWidth', LINE_WIDTH)
hold off

xlabel('ms', 'FontSize', FONT_SIZE)
ylabel('Amplitude (uV)', 'FontSize', FONT_SIZE)
title('S001, tip-h, 80 dB, n=500 single traces', 'FontSize', FONT_SIZE)