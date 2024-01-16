function plotstack_average_julienCAPABR(varargin)
% Plot a stack of average CAPs and ABRs using single traces, one stack per frequency. Saves
% plots in multiple image formats.
%
% Default is for user to select file in dialog box when no input parameters
% are specified.
%
% Optional input syntax for specifying input file: plotstack_averageABR(path, filename)
%
% TODO: 
% [ ] Some TDT systems save average ABR CSV file in units of V instead of
% nV. Make import_averageABRcsv.m check units of ABR voltages and convert
% values to nV if necessary. (Added 10/6/21).
%
% 7/12/2023 George Liu
% Dependencies: import_ABRcsv.m, merge_singletraceABR_polarities, savefig_multiformat.m

%% Constants
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures
% SAMPLE_PERIOD = 40.96; % sample period of 40.96 us/sample in average trace
SAMPLE_PERIOD = 5.12; % sample period of 5.12 us/sample in single trace
SAMPLE_PERIOD_MS = SAMPLE_PERIOD/1000; % sample period in ms

NUM_CHANNELS = 2; % 1 = ABR, 2 = CAP
CHANNEL_KEY = {'ABR', 'CAP'};
CHANNEL_YLIM = {[-1000, 1000], [-1000, 1000]};

NUM_LEVELS = 9;
% N_FREQ = 2; % 16 kHz, 32 kHz

PLOT_FIGURES = 1; % Hard code logical for whether or not to plot figures 
SAVE_FIGURES = 0; % Hard code logical for whether or not to save figures 

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

% Analyze metadata in single trial ABR/CAP filename
% Filename example: 20230628_tmie2_ABRCAP_pre-0-55-2-1.csv
% group 0
% SGI 55
% channel 2
% run (?) 1
dash_loc = strfind(filename, '-');
dot_loc = strfind(filename, '.');

num_dashes = length(dash_loc);
group_ind = dash_loc(num_dashes - 3) + 1 : dash_loc(num_dashes - 2) - 1;
SGI_ind = dash_loc(num_dashes - 2) + 1 : dash_loc(num_dashes - 1) - 1;
channel_ind = dash_loc(num_dashes - 1) + 1 : dash_loc(num_dashes) - 1;
run_ind = dash_loc(num_dashes) + 1 : dot_loc - 1;

% Metadata of selected file
group = str2double(filename(group_ind));
% SGI = str2double(filename(SGI_ind));
% channel = str2double(filename(channel_ind)); % 1 = ABR, 2 = CAP
run = str2double(filename(run_ind));

% Obtain all single trial ABR/CAP files from same "run" / timepoint, which
% have the same filename stem.
filename_stem = filename(1 : dash_loc(end-3));
listing = dir(fullfile(path, [filename_stem, '*.csv']));

allnames = {listing.name}'; % column cell array of filenames

extract = @(C, k) cellfun(@(c) str2double((c(k))), C) ; % cell function to extract k'th element from each cell in cell array C.
all_groups = extract(allnames, group_ind);
all_SGI = extract(allnames, SGI_ind);
all_channels = extract(allnames, channel_ind); % 1 = ABR, 2 = CAP
all_runs = extract(allnames, run_ind);

for i = 1:NUM_CHANNELS
    disp(['Working on channel ', num2str(i), ' out of ', num2str(NUM_CHANNELS), ': ', CHANNEL_KEY{i}])
    % Obtain ABR/CAP traces, ordered by SGI (same group and run as selected trace)
    is_data = all_channels == i;
    is_group = all_groups == group;
    is_run = all_runs == run;

    data_filenames = allnames(is_data & is_group & is_run);
    [~, sgi_order] = sort(all_SGI(is_data & is_group & is_run));
    data_filenames_ordered = data_filenames(sgi_order);
    
    num_data = length(data_filenames_ordered);
%     num_levels = num_data / N_FREQ; % assume only 2 tone pip frequencies, 16 and 32 kHz
    n_freq = num_data / NUM_LEVELS; % assume only 10 stimulus levels, 0, 10, ..., 90 dB
    
    if PLOT_FIGURES
        fig = figure;
        t = tiledlayout(NUM_LEVELS, n_freq, 'TileIndexing', 'columnmajor');
    end
        
    % Collect all average ABR traces from single trial ABR recordings
    y_avg_cache = cell(num_data, 1);
    A_csv_cache = cell(num_data, 1);
    freq_csv_cache = cell(num_data, 1);
    for j = 1:num_data
        this_filename = data_filenames_ordered{j};
        disp(['Working on file ', num2str(j), ' out of ', num2str(num_data), ': ', this_filename])
        [X_csv, A_csv, freq_csv] = import_ABRcsv(this_filename, path);
    
        % merge single trace pairs with alternating polarity to cancel cochlear
        % microphonic heterogeneity
        y = merge_singletraceABR_polarities(X_csv); % 1026 -> 513 single traces. Matrix of size (SAMPLES, m_traces/2)
        y_avg = mean(y, 2);
        
        % Cache variables
        y_avg_cache{j, 1} = y_avg;
        A_csv_cache{j, 1} = A_csv;
        freq_csv_cache{j, 1} = freq_csv;
        
        num_samples = size(y, 1);
%         n_traces = size(y, 2);
        X = SAMPLE_PERIOD_MS * (1:num_samples); % time in ms
    end
        
    % Make sure bounds of plot are within y limits
    ylim_max = CHANNEL_YLIM{i};
    ymax_data = max(cell2mat(y_avg_cache));
    ymin_data = min(cell2mat(y_avg_cache));
    if ymax_data > ylim_max(2)
        ylim_max(2) = ymax_data;
    end
    if ymin_data < ylim_max(1)
        ylim_max(1) = ymin_data;
    end
            
    for j = 1:num_data
        if PLOT_FIGURES
            y_avg = y_avg_cache{j};
            A_csv = A_csv_cache{j};
            freq_csv = freq_csv_cache{j};
            
            nexttile(t)
            plot(X, y_avg, 'LineWidth', 3)
            ylim(ylim_max)
                        
            % show ylabels for first column only
            if j/NUM_LEVELS <= 1
                ylabel([num2str(A_csv), ' dB (nV)'], 'FontSize', 24)
                % rotate y label to make horizontal
                hYLabel = get(gca,'YLabel');
                set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right', 'FontSize', 24)
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

            % Remove y-axis labels from all but first column
            if j/NUM_LEVELS > 1
                set(gca,'yticklabel', [])
            end

            % Set text size of current axis
            set(gca,'FontSize',24)

            % Show frequency as title above top-most tile of each column
            if mod(j, NUM_LEVELS)==1
                title([num2str(freq_csv), ' Hz']) 
            end
        end
    
    
    end
    
    xlabel(t, 'Time (ms)', 'FontSize', 24)
%     ylabel(t, 'Voltage (nV)', 'FontSize', 24)
    title(t, CHANNEL_KEY{i}, 'FontSize', 24)
    
    t.TileSpacing = 'none';
    t.Padding = 'tight';

    % Maximize figure window size
    fig.WindowState = 'maximized';

    if SAVE_FIGURES
        % Save figure
        %     disp('Saving figure')
        [~, save_file, ~] = fileparts(filename);
        savefig_multiformat(gcf, SAVE_PATH, [save_file, '_', CHANNEL_KEY{i}])
    end

end

% % Obtain CAP traces, ordered by SGI (same group and run as selected trace)
% is_cap = all_channels == CAP_CHANNEL;
% 
% cap_filenames = allnames(is_cap & is_group & is_run);
% [~, sgi_order2] = sort(all_SGI(is_cap & is_group & is_run));
% cap_filenames_ordered = cap_filenames(sgi_order2);

% [M, A_csv, freq_csv] = import_averageABRcsv(filename, path);
% num_samples = size(M, 2);
% X = SAMPLE_PERIOD_MS * (1:num_samples); % time in ms
% A_levels = unique(A_csv);
% num_levels = length(A_levels);

disp('Done')

end