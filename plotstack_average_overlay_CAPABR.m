% function plotstack_average_overlay_CAPABR(varargin)
% Plot average CAP and ABR data with overlay of pre-injection traces. 
%
% 8/4/2023 George Liu
% Dependencies: none

%% Constants
LOAD_PATH = 'd:\users\admin\Documents\George\Results';
SAVE_PATH = 'd:\users\admin\Documents\George\Results_analyzed'; % path for saving figures

SAMPLE_PERIOD = 5.12; % sample period of 5.12 us/sample in single trace
SAMPLE_PERIOD_MS = SAMPLE_PERIOD/1000; % sample period in ms

SAVE_FIGURES = 0;
NUM_CHANNELS = 2; % 1 = ABR, 2 = CAP
CHANNEL_KEY = {'ABR', 'CAP'};

CIRCLE_SIZE = 54; % default 36

%% Load data
disp('Select file of post injection MAT results')
[filename,path] = uigetfile('*.mat', 'Select post injection MAT file');
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

date = substrings{1};
label = substrings{2}{1};
timepoint = substrings{3}{1};
% ABRCAP = substring{4}{1};

if strcmp(timepoint, 'pre')
    disp(' Skipping pre timepoint')
    return
end

filename_prestem = [date, '_', label, '_'];
my_vars = {'y_avg_cache_descending', 'y_avg_cm_cache_descending', ...
    'all_freqs', 'A_descending', 'threshold', 'metric'};

for i = 1:NUM_CHANNELS
    filename_poststem = ['_', CHANNEL_KEY{i}, '_results.mat'];
    filename_pre = [filename_prestem, 'pre', filename_poststem];
%     filename_pre = [filename_prestem, 'post10', filename_poststem];
    filename_post = [filename_prestem, timepoint, filename_poststem];
    
    % Load selected and baseline ABR / CAP average trace data
    load(fullfile(LOAD_PATH, filename_post), my_vars{:})
    y_avg_sel = y_avg_cache_descending; % each cell entry is sorted stack of average traces for a frequency
    y_cm_sel = y_avg_cm_cache_descending;

%     filename_pre = [filename_prestem, 'pre', filename_poststem]; % REMOVE
    
    load(fullfile(path, filename_pre), my_vars{:})
    y_avg_pre = y_avg_cache_descending;
    y_cm_pre = y_avg_cm_cache_descending;
    
    % Obtain plotting parameters
    n_freq = length(all_freqs);
    n_levels = length(A_descending);
    num_samples = size(y_avg_sel{1}, 2);
    X = SAMPLE_PERIOD_MS * (1:num_samples); % time in ms
    
    % Make sure bounds of plot are within y limits
    ylim_max = [];
    ymax_data = max(cell2mat(y_avg_sel), [], "all");
    ymin_data = min(cell2mat(y_avg_sel), [], "all");
    ymax_data_pre = max(cell2mat(y_avg_pre), [], "all");
    ymin_data_pre = min(cell2mat(y_avg_pre), [], "all");
	ylim_max(2) = max([ymax_data, ymax_data_pre]);
    ylim_max(1) = min([ymin_data, ymin_data_pre]);
    
    % Plot
    fig = figure;
    t = tiledlayout(n_levels, n_freq, 'TileIndexing', 'columnmajor');
    
    for j = 1:n_freq
        for k = 1:n_levels
            nexttile(t)
%             plot(X, y_avg_pre{j}(k, :), 'LineWidth', 2.5, 'Color', [0 0.4470 0.7410])
            hold on
            plot(X, y_avg_sel{j}(k, :), 'LineWidth', 2.5, 'Color', [0.8500 0.3250 0.0980])
            
            % Calculate wave 1 amplitude and latency 
            y_avg = y_avg_sel{j}(k, :);
            if k==1
                % Get wave 1 peak and following trough at highest stimulus
                % level
                if i==1 % ABR
                    [peak_pt, trough_pt, amp, lat_peak, lat_trough] = get_wave1_averageABR(X, y_avg);
                else % CAP
                    [peak_pt, trough_pt, amp, lat_peak, lat_trough] = get_wave1_averageCAP(X, y_avg);
                end
            else
                % Get wave 1 peak and following trough at lower stimulus
                % level. Ensure peak latency does not decrease.
                if i==1 % ABR
                    [peak_pt, trough_pt, amp, lat_peak, lat_trough] = get_wave1_averageABR(X, y_avg, peak_pt(1), trough_pt(1));
                else % CAP
                    [peak_pt, trough_pt, amp, lat_peak, lat_trough] = get_wave1_averageCAP(X, y_avg, peak_pt(1), trough_pt(1));
                end
            end
            
            scatter(peak_pt(1), peak_pt(2), CIRCLE_SIZE, 'green', 'filled') % peak
            scatter(trough_pt(1), trough_pt(2), CIRCLE_SIZE, 'red', 'filled') % trough
            
%             if j==1 && k==1 && i==2
%                 % Calculate CAP wave 1 amplitude
%                 [peak_pt, trough_pt, amp, latency, latency_trough] = get_wave1_averageABR(X, -1*y_avg_pre{1}(1, :));
%                 [peak_pt2, trough_pt2, amp2, latency2, latency_trough2] = get_wave1_averageABR(X, -1*y_avg_sel{1}(1, :), peak_pt(1), trough_pt(1));
% %                 [peak_pt2, trough_pt2, amp2, latency2, latency_trough2] = get_wave1_averageABR(X, -1*y_avg_sel{1}(1, :), 2.1, 2.7);
%                 
%                 % Show CAP wave 1 peak and trough
%                 scatter(peak_pt(1), -1*peak_pt(2), CIRCLE_SIZE, 'green', 'filled') % peak
%                 scatter(trough_pt(1), -1*trough_pt(2), CIRCLE_SIZE, 'red', 'filled') % trough
%                 cap_amp_pre = abs(trough_pt(2) - peak_pt(2));
%                 scatter(peak_pt2(1), -1*peak_pt2(2), CIRCLE_SIZE, 'green') % peak
%                 scatter(trough_pt2(1), -1*trough_pt2(2), CIRCLE_SIZE, 'red') % trough
%                 cap_amp_post = abs(trough_pt2(2) - peak_pt2(2));
% %                 disp(['Pre CAP wave 1 amp (nv): ', num2str(cap_amp_pre)])
%                 disp(['Post 10 CAP wave 1 amp (nv): ', num2str(cap_amp_pre)])
%                 disp(['Post 60 CAP wave 1 amp (nv): ', num2str(cap_amp_post)])
%                 disp(['Change CAP wave 1 amp (nv): ', num2str(cap_amp_post - cap_amp_pre)])
%             end
            
            hold off
            ylim(ylim_max)
                        
            % show ylabels for first column only
            if j == 1
                ylabel([num2str(A_descending(k)), ' dB (nV)'], 'FontSize', 24)
                % rotate y label to make horizontal
                hYLabel = get(gca,'YLabel');
                set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right', 'FontSize', 24)
            end

            % Adjust plot appearance
            % Remove extraneous axis tick marks and x-axis from all but bottom
            % tile
            set(gca,'box','off')
            if mod(k, n_levels) ~= 0 % 
                set(gca,'xtick',[])
                set(gca,'xticklabel',[])
                set(gca,'XColor','none')
            end

            % Remove y-axis labels from all but first column
            if j > 1
                set(gca,'yticklabel', [])
            end

            % Set text size of current axis
            set(gca,'FontSize',24)

            % Show frequency as title above top-most tile of each column
            if k==1
                title([num2str(all_freqs(j)), ' Hz']) 
            end
        end
    end
    
    % Format figure
    xlabel(t, 'Time (ms)', 'FontSize', 24)
%     ylabel(t, 'Voltage (nV)', 'FontSize', 24)
    title(t, CHANNEL_KEY{i}, 'FontSize', 24)

    t.TileSpacing = 'none';
    t.Padding = 'tight';
    
    % Legend
    lg = legend({'Pre', timepoint}, 'Orientation', 'Vertical');
    lg.Layout.Tile = 'East'; % <-- Legend placement with tiled layout

    % Maximize figure window size
    fig.WindowState = 'maximized';    

    % Save figure
    if SAVE_FIGURES 
        save_file = [filename_prestem, timepoint];
        savefig_multiformat(gcf, SAVE_PATH, [save_file, '_overlaypre_', CHANNEL_KEY{i}])
    end
    
    %% Calculate CM if CAP
    if i==2
        % Calculate CM
        max_8khz_cm_post  = y_cm_sel{1}(1, 1:floor(num_samples/2));
        max_8khz_cm_pre  = y_cm_pre{1}(1, 1:floor(num_samples/2));
        amp_cm_post = max(max_8khz_cm_post) - min(max_8khz_cm_post);
        amp_cm_pre = max(max_8khz_cm_pre) - min(max_8khz_cm_pre);
        amp_cm_change = amp_cm_post - amp_cm_pre;
%         disp(['Pre cochlear microphonic amp (nV): ', num2str(amp_cm_pre)])
        disp(['Post 10 cochlear microphonic amp (nV): ', num2str(amp_cm_pre)])
        disp(['Post 60 cochlear microphonic amp (nV): ', num2str(amp_cm_post)])
        disp(['Change cochlear microphonic amp (nV): ', num2str(amp_cm_change)])
        

    end
end
    
    