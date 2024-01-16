% function pca_singletraces_CAPABR(varargin)
% Plot and determine PCA of single traces file of CAP or ABR at one stimulus level and frequency. Saves
% plots in multiple image formats.
%
% Default is user to select file in dialog box when no input parameters
% are specified.
%
% 7/11/2023 George Liu
% Dependencies: import_ABRcsv.m, merge_singletraceABR_polarities,
% savefig_multiformat.m, get_fft_CAPABR.m

opengl('save', 'software') % prevent crashing due to low-level graphics error using Sarah Office computer's graphics card

%% Constants
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures
% SAMPLE_PERIOD = 40.96; % sample period of 40.96 us/sample in average trace
SAMPLE_PERIOD = 5.12; % sample period of 5.12 us/sample in single trace
SAMPLE_PERIOD_MS = SAMPLE_PERIOD/1000; % sample period in ms

CHANNEL_KEY = {'ABR', 'CAP'};
% YLIM_DEFAULT = [-1000, 1000];

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

num_dashes = length(dash_loc);
channel_ind = dash_loc(num_dashes - 1) + 1 : dash_loc(num_dashes) - 1;
channel = str2double(filename(channel_ind));

% Load data
[X_csv, A_csv, freq_csv] = import_ABRcsv(filename, path);
    
% merge single trace pairs with alternating polarity to cancel cochlear
% microphonic heterogeneity
y = merge_singletraceABR_polarities(X_csv); % 1026 -> 513 single traces. Matrix of size (SAMPLES, m_traces/2)
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

fig = figure;

% Overlay single traces
n_traces = size(y, 2);
cmap = spring(n_traces);
hold on
for i = 1:n_traces
    plot(X, X_csv(:, i), 'Color', [cmap(i, 1), cmap(i, 2), cmap(i, 3), 0.3]) % transparency set to 4th color value
%     plot(X, X_csv(:, i), 'Color', [0, 0, 0, 0.2]) % transparency set to 4th color value
end

plot(X, y_avg, 'blue', 'LineWidth', 3)
% ylim(ylim_max)
hold off

ylabel([num2str(A_csv), ' dB (nV)'], 'FontSize', 24)
% rotate y label to make horizontal
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right', 'FontSize', 24)

% Adjust plot appearance
set(gca,'box','off')

% Set text size of current axis
set(gca,'FontSize',24)

% Show frequency as title above top-most tile of each column
title([CHANNEL_KEY{channel}, ' @ ', num2str(freq_csv), ' Hz']) 
xlabel('Time (ms)', 'FontSize', 24)

% Maximize figure window size
fig.WindowState = 'maximized';

if SAVE_FIGURES
    % Save figure
    %     disp('Saving figure')
    [~, save_file, ~] = fileparts(filename);
    savefig_multiformat(gcf, SAVE_PATH, [save_file, '_', num2str(A_csv), 'dB_', num2str(freq_csv/1000), 'kHz_', CHANNEL_KEY{channel}, '_singletraces'])
end


%% Plot representative sample of single traces
NUM_COLS = 5;
NUM_ROWS = 5;
fig2 = figure;
t = tiledlayout(NUM_ROWS , NUM_COLS, 'TileIndexing', 'rowmajor');

trace_ind = randi(n_traces, NUM_ROWS, NUM_COLS);
trace_ind = reshape(sort(trace_ind(:)), NUM_ROWS, NUM_COLS);
for i = 1:NUM_ROWS
    for j = 1:NUM_COLS
        nexttile(t)
        plot(X, X_csv(:, trace_ind(i, j)), 'blue', 'LineWidth', 3) 
        
        title(['Trace ', num2str(trace_ind(i, j))], 'FontSize', 24) 
    end
end
title(t, [CHANNEL_KEY{channel}, ' @ ', num2str(freq_csv), ' Hz'], 'FontSize', 24)
xlabel(t, 'Time (ms)', 'FontSize', 24)
ylabel(t, [num2str(A_csv), ' dB (nV)'], 'FontSize', 24, 'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right')

fig2.WindowState = 'maximized'; % Maximize figure window size

if SAVE_FIGURES
    % Save figure
    %     disp('Saving figure')
    [~, save_file, ~] = fileparts(filename);
    savefig_multiformat(gcf, SAVE_PATH, [save_file, '_', num2str(A_csv), 'dB_', num2str(freq_csv/1000), 'kHz_', CHANNEL_KEY{channel}, '_singletracesexamples'])
end

%% Calculate PCA of single ABR traces at single frequency and stimulus 
% level for noise reduction and dimension reduction
q = 7;
[coeff, score, latent, tsquared, explained, mu] = pca(y', 'NumComponents', q);

% Compute and display explained variation
text_topcompsvar = strcat("Top ", string(q), " principle components explain ", ...
    sprintf('%.3g', sum(explained(1:q))), "% of variation");
disp(text_topcompsvar)

c = parula(q);
fig3 = figure;
legend_labels = cell(q + 1, 1);
hold on
for i = 1:q 
    plot(X, coeff(:, i), 'Color', c(i, :), 'LineWidth', 3)
    legend_labels{i, 1} = ['Comp ', num2str(i), ', var=', num2str(explained(i))];
end
% Plot average trace
max_coeff = max(coeff(:, 1:q), [], 'all');
y_avg_norm = y_avg/max(y_avg)*max_coeff;
plot(X, y_avg_norm, 'k--', 'LineWidth', 3)
legend_labels{q + 1, 1} = 'Normalized average';
hold off

% Set text size of current axis
set(gca,'FontSize',24)

legend(legend_labels)
title(strcat(text_topcompsvar, [', ', CHANNEL_KEY{channel}, ' @ ', num2str(freq_csv), ' Hz']), 'FontSize', 24)
xlabel('Time (ms)', 'FontSize', 24)
ylabel([num2str(A_csv), ' dB (nV)'], 'FontSize', 24, 'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right')

fig3.WindowState = 'maximized'; % Maximize figure window size

if SAVE_FIGURES
    % Save figure
    %     disp('Saving figure')
    [~, save_file, ~] = fileparts(filename);
    savefig_multiformat(gcf, SAVE_PATH, [save_file, '_', num2str(A_csv), 'dB_', num2str(freq_csv/1000), 'kHz_', CHANNEL_KEY{channel}, '_pca', num2str(q)])
end


%% Plot FFT of each PCA component
fig4 = figure;
t2 = tiledlayout(q, 2, 'TileIndexing', 'rowmajor');
for i = 1:q
    nexttile(t2)
    plot(X, coeff(:, i), 'Color', c(i, :), 'LineWidth', 3)
    hold on
    plot(X, y_avg_norm, 'k--', 'LineWidth', 3)
    hold off
    set(gca,'FontSize',24)
    ylabel([num2str(i), '(', sprintf('%.3g', explained(i)), '%)'], 'FontSize', 24, 'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right')
    
    % Show xlabel for bottom-most tile of each column
    if i==q
        xlabel('Time (ms)', 'FontSize', 24)
    end
    
    % Show title for top-most tile of each column
    if i==1
        title('Component', 'FontSize', 24)
    end
    
    nexttile(t2)
    [f, P] = get_fft_CAPABR(coeff(:, i));
    plot(f, P, 'Color', c(i, :), 'LineWidth', 3)
    xlim([0, 3000])
    set(gca,'FontSize',24)
    
    % Show xlabel for bottom-most tile of each column
    if i==q
        xlabel("f (Hz)")
    end
    
    % Show title for top-most tile of each column
    if i==1
        title('FFT', 'FontSize', 24)
    end
    
    % Remove y-axis labels from all but first column
    set(gca,'yticklabel', [])

end

title(t2, ['PCA of ', CHANNEL_KEY{channel}, ' @ ', num2str(freq_csv), ' Hz, ',  num2str(A_csv), ' dB'], 'FontSize', 24)
ylabel(t2, 'nV', 'FontSize', 24)
t2.TileSpacing = 'none';
t2.Padding = 'tight';
fig4.WindowState = 'maximized'; % Maximize figure window size

if SAVE_FIGURES
    % Save figure
    %     disp('Saving figure')
    [~, save_file, ~] = fileparts(filename);
    savefig_multiformat(gcf, SAVE_PATH, [save_file, '_', num2str(A_csv), 'dB_', num2str(freq_csv/1000), 'kHz_', CHANNEL_KEY{channel}, '_pca', num2str(q), '_fft'])
end

disp('Done')

% end