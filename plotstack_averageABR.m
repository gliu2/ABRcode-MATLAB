function outTable = plotstack_averageABR(varargin)
% Plot a stack of average ABRs in CSV file, one stack per frequency. Saves
% plots in multiple image formats.
% Analyzes wave 1 amplitude and latency in highest stimulus waveform.
% Estimate ABR threshold.
%
% Optional input: plotstack_averageABR(path, filename)
%
% 7/11/2021 George Liu
% Dependencies: import_averageABRcsv.m, savefig_multiformat.m, get_wave1_averageABR

%% Constants
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures
SAMPLE_PERIOD = 40.96; % sample period of 40.96 us/sample in average trace
TIME_WINDOW = [1, 2]; % time in ms to analyze portion of ABR trace with inner product method

%% Load average trace data
if nargin == 2
    path = varargin{1};
    filename = varargin{2};
elseif nargin == 0
    [filename,path] = uigetfile('*.csv');
else
    disp('Warning: number of input arguments is not 0 or 2!')
    [filename,path] = uigetfile('*.csv');
end

disp(['Opening ', fullfile(path, filename)])
[M, A_csv, freq_csv] = import_averageABRcsv(filename, path);
num_samples = size(M, 2);
X = SAMPLE_PERIOD / 1000 * (1:num_samples); % time in ms
A_levels = unique(A_csv);
num_levels = length(A_levels);

% Analyze ABRs for each frequency separately
[freq_unique, ~, ic] = unique(freq_csv);
n_freq = length(freq_unique);

% Initialize cache variables
Frequency = freq_unique;
Stimulus_amplitude = zeros(n_freq, num_levels);
Threshold_liberman = zeros(n_freq, 1);
Threshold_oghalai = zeros(n_freq, 1);
Threshold_innerprod = zeros(n_freq, 1);
Threshold_wave1amp = zeros(n_freq, 1);
Metric_liberman = zeros(n_freq, num_levels);
Metric_oghalai = zeros(n_freq, num_levels);
Metric_innerprod = zeros(n_freq, num_levels);
Metric_innerprod_auc = zeros(n_freq, num_levels);
Metric_wave1amp = zeros(n_freq, num_levels);
Metric_wave1lat = zeros(n_freq, num_levels);

% Obtain y axis scale
ylim_max = [-1200, 1200];
% Check to make sure bounds of plot are within y limits
ymax_data = max(M(:));
ymin_data = min(M(:));
if ymax_data > ylim_max(2)
    ylim_max(2) = ymax_data;
end
if ymin_data < ylim_max(1)
    ylim_max(1) = ymin_data;
end

% Plot ABRs in column stacks for each frequencies, in one tiled layout
fig = figure;
t = tiledlayout(num_levels, n_freq, 'TileIndexing', 'columnmajor');
for ff = 1:n_freq
    this_freq = freq_unique(ff);
%     disp(['Working on ', num2str(this_freq), ' Hz...'])
    
    % amplitudes and DPOAE values for this frequency only
    A_csv2 = A_csv(ic==ff, :);
    M2 = M(ic==ff, :);
    
%     num_levels = size(M2, 1);
    dB = cell(num_levels, 1);
    dB(:) = {'dB (nV)'};

    [A_descending, I] = sort(A_csv2, 'descend');
    M_sorted = M2(I, :);
    A_descending_cell = cellstr(num2str(A_descending));
    newYlabels=cellfun(@(x,y) [x  ' ' y], A_descending_cell, dB, 'un', 0);

    %% Plot ABRs in vertical stack
%     figure
%     t = tiledlayout(num_levels, 1);
%     s = stackedplot(X, M_sorted', 'Title', plot_title, 'DisplayLabels', newYlabels, 'LineWidth', 1.5);
    
    % Calculate threshold
    [thresh_cache_liberman, this_metric_liberman] = get_thresh_averageABR_liberman(M_sorted, A_descending);
    [thresh_cache_oghalai, this_metric_oghalai] = get_thresh_averageABR_oghalai(M_sorted, A_descending, false);

    % initialize variables
    is_abovenoise = zeros(num_levels, 1);
    for i = 1:num_levels
        nexttile(t)
        y = M_sorted(i, :);
        plot(X, y, 'LineWidth', 1.5)
        ylim(ylim_max)
        
        % show ylabels for first column only
        if ff==1
            ylabel(newYlabels{i})
            % rotate y label to make horizontal
            hYLabel = get(gca,'YLabel');
            set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right')
        end

        % Mark wave 1 peak and following trough
        if i==1
            % Get wave 1 peak and following trough at highest stimulus
            % level
            [peak_pt, trough_pt, amp, lat] = get_wave1_averageABR(X, y);
        else
            % Get wave 1 peak and following trough at lower stimulus
            % level. Ensure peak latency does not decrease.
            [peak_pt, trough_pt, amp, lat] = get_wave1_averageABR(X, y, peak_pt(1), trough_pt(1));
        end
        hold on
        scatter(peak_pt(1), peak_pt(2), [], 'green') % peak
        scatter(trough_pt(1), trough_pt(2), [], 'red') % trough
        hold off
        
        % Display wave 1 measurements in corner of plot
        xL=xlim;
        yL=ylim;
        str = sprintf('P-P_I = %.0f nV', amp);
        if A_descending(i) > thresh_cache_liberman
            text(0.99*xL(2), 0.99*yL(1), str,'HorizontalAlignment','right','VerticalAlignment','bottom', 'FontWeight', 'bold')
            is_abovenoise(i,1)=1;
        else
            text(0.99*xL(2), 0.99*yL(1), str,'HorizontalAlignment','right','VerticalAlignment','bottom')
            is_abovenoise(i,1)=0;
        end
        
        % Remove extraneous axis tick marks and x-axis from all but bottom
        % tile
        set(gca,'box','off')
        if i~=num_levels
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
        end
        
        % Remove y-axis labels from all but first column
        if ff~=1
            set(gca,'yticklabel',[])
        end
        
        % Set text size of current axis
        set(gca,'FontSize',10)
        
        % Show frequency as title above top-most tile of each column
        if i==1
            title([num2str(this_freq), ' Hz'])
        end
        
        % Cache wave 1 amplitude and latency for highest stimulus level
        Metric_wave1amp(ff, i) = amp; % nV
        Metric_wave1lat(ff, i) = lat; % ms
    end

    xlabel('Time (ms)')
    
    % Calculate threshold based on wave 1 amplitude
    noise_std = std(M_sorted(end, :));
    wave1_thresh = get_thresh_wave1amp_oghalai(Metric_wave1amp(ff, :), A_descending, noise_std, false);
    
    % Cache variables
    Threshold_wave1amp(ff) = wave1_thresh;
    Stimulus_amplitude(ff, :) = A_descending';
    Threshold_liberman(ff) = thresh_cache_liberman; % dB
    Threshold_oghalai(ff) = thresh_cache_oghalai;
    Metric_liberman(ff, :) = this_metric_liberman;
    Metric_oghalai(ff, :) = this_metric_oghalai;
end

t.TileSpacing = 'none';
t.Padding = 'tight';
%     plot_title = ['ABR @ ', num2str(this_freq), ' Hz'];
%     title(t, plot_title)
%     ylabel(t, 'Amplitude (nV)')

% Maximize figure window size
fig.WindowState = 'maximized';

% % Save figure
% %     disp('Saving figure')
% [~, save_file, ~] = fileparts(filename);
% savefig_multiformat(gcf, SAVE_PATH, save_file)

% Obtain single trace inner product threshold and AUC level functions
[auc_cache, ip_mean_cache, ip_ste_cache, thresh_cache_innerprod, auc_lower_cache, auc_upper_cache] = get_roc_innerprod_ABR(fullfile(path, filename));
Threshold_innerprod = thresh_cache_innerprod;
Metric_innerprod = ip_mean_cache;
Metric_innerprod_ste = ip_ste_cache;
Metric_innerprod_auc = auc_cache;
Metric_innerprod_auc_lower = auc_lower_cache;
Metric_innerprod_auc_upper = auc_upper_cache;

% Obtain wave 1 time window (1-2 ms) single trace inner product threshold and AUC level functions
[auc_window_cache, ip_window_mean_cache, ip_window_ste_cache, thresh_window_cache_innerprod, auc_window_lower_cache, auc_window_upper_cache] = get_roc_innerprod_ABR(fullfile(path, filename), false, TIME_WINDOW);
Threshold_innerprod_window = thresh_window_cache_innerprod;
Metric_innerprod_window = ip_window_mean_cache;
Metric_innerprod_window_ste = ip_window_ste_cache;
Metric_innerprod_auc_window = auc_window_cache;
Metric_innerprod_auc_window_lower = auc_window_lower_cache;
Metric_innerprod_auc_window_upper = auc_window_upper_cache;

% Write cached variables to excel table
% Structure of table is: filename, frequency, method, threshold, metric
% level function.
% Listing of filenames
[~, baseFileNameNoExt, ~] = fileparts(filename);
Filenames = cell(n_freq, 1);
Filenames(:) = {baseFileNameNoExt};  
Filenames = convertCharsToStrings(Filenames);

% Metric names
Metric_names_amplitude = repmat("Amplitude", n_freq, 1);
Metric_names_liberman = repmat("Liberman", n_freq, 1);
Metric_names_oghalai = repmat("Oghalai", n_freq, 1);
Metric_names_innerprod = repmat("Innerprod", n_freq, 1);
Metric_names_innerprod_auc = repmat("Innerprod_auc", n_freq, 1);
Metric_names_wave1amp = repmat("Wave1amp", n_freq, 1);
Metric_names_wave1lat = repmat("Wave1lat", n_freq, 1);
blank_space = nan(n_freq, 1);

% Stimulus amplitudes in place of metrics
outTable_amplitudes = table(Filenames, Frequency, Metric_names_amplitude, blank_space , ...
    Stimulus_amplitude);
outTable_amplitudes.Properties.VariableNames{'Metric_names_amplitude'} = 'Metric';
outTable_amplitudes.Properties.VariableNames{'blank_space'} = 'Threshold';
outTable_amplitudes.Properties.VariableNames{'Stimulus_amplitude'} = 'Metric_values';
% Liberman method results
outTable_liberman = table(Filenames, Frequency, Metric_names_liberman, Threshold_liberman, ...
    Metric_liberman);
outTable_liberman.Properties.VariableNames{'Metric_names_liberman'} = 'Metric';
outTable_liberman.Properties.VariableNames{'Threshold_liberman'} = 'Threshold';
outTable_liberman.Properties.VariableNames{'Metric_liberman'} = 'Metric_values';
% Oghalai method results
outTable_oghalai = table(Filenames, Frequency, Metric_names_oghalai, Threshold_oghalai, ...
    Metric_oghalai);
outTable_oghalai.Properties.VariableNames{'Metric_names_oghalai'} = 'Metric';
outTable_oghalai.Properties.VariableNames{'Threshold_oghalai'} = 'Threshold';
outTable_oghalai.Properties.VariableNames{'Metric_oghalai'} = 'Metric_values';
% Inner product method results 
outTable_innerprod = table(Filenames, Frequency, Metric_names_innerprod, Threshold_innerprod, ...
    Metric_innerprod);
outTable_innerprod.Properties.VariableNames{'Metric_names_innerprod'} = 'Metric';
outTable_innerprod.Properties.VariableNames{'Threshold_innerprod'} = 'Threshold';
outTable_innerprod.Properties.VariableNames{'Metric_innerprod'} = 'Metric_values';
% Inner product AUC method results
outTable_innerprod_AUC = table(Filenames, Frequency, Metric_names_innerprod_auc, Threshold_innerprod, ...
    Metric_innerprod_auc);
outTable_innerprod_AUC.Properties.VariableNames{'Metric_names_innerprod_auc'} = 'Metric';
outTable_innerprod_AUC.Properties.VariableNames{'Threshold_innerprod'} = 'Threshold';
outTable_innerprod_AUC.Properties.VariableNames{'Metric_innerprod_auc'} = 'Metric_values';
% Wave 1 amplitude method results
outTable_wave1amp = table(Filenames, Frequency, Metric_names_wave1amp, Threshold_wave1amp, ...
    Metric_wave1amp);
outTable_wave1amp.Properties.VariableNames{'Metric_names_wave1amp'} = 'Metric';
outTable_wave1amp.Properties.VariableNames{'Threshold_wave1amp'} = 'Threshold';
outTable_wave1amp.Properties.VariableNames{'Metric_wave1amp'} = 'Metric_values';
% Wave 1 latency
outTable_wave1lat = table(Filenames, Frequency, Metric_names_wave1lat, Threshold_wave1amp, ...
    Metric_wave1lat);
outTable_wave1lat.Properties.VariableNames{'Metric_names_wave1lat'} = 'Metric';
outTable_wave1lat.Properties.VariableNames{'Threshold_wave1amp'} = 'Threshold';
outTable_wave1lat.Properties.VariableNames{'Metric_wave1lat'} = 'Metric_values';

outTable = [outTable_amplitudes; outTable_liberman; outTable_oghalai; outTable_innerprod; outTable_innerprod_AUC; ...
    outTable_wave1amp; outTable_wave1lat];

% % Save excel sheet output
% table_filename = [baseFileNameNoExt, '_output.xlsx'];
% writetable(outTable, fullfile(SAVE_PATH, table_filename), 'Sheet', 1, 'Range', 'A1')

% disp('Done')

end