% Wrapper script to pre, post24hour, and post1week ABR files for one mouse.
% Quickly iterate testing of single trial ABR analysis method.
%
% SAVE_PATH automatically set to "Results" folder to save output files
%
% Dependencies: plotstack_averageABR.m
%
% 12/8/21 George Liu
% Last edit: 9/20/2022

close all
opengl('save', 'software') % prevent crashing due to low-level graphics error using Sarah Office computer's graphics card
addpath 'D:\users\admin\Documents\GitHub\ABRcode-MATLAB\RAP__Risk_Assessment_Plot_' %path to DeLong AUC comparison code

%% Constants
IS_ABR = 1; % MANUALLY SET TO 1 (ABR) OR 0 (DPOAE)
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures
LOAD_PATH = 'd:\users\admin\Documents\George\finalTable_metadata_4-18-22\finalTable_metadata_9-13-22.mat'; % load variable finalTable_metadata
% ABR_PATH = '\\Ricci-abr\d\George\ABR'; % path to data saved on ABR computer - slow to load in data
ABR_PATH = 'D:\George-abr\ABR'; % path to local copy of data, 12-28-21

% HARD CODE WHICH MOUSE AND DATES TO ANALYZE
% Example mice:
%   - 'b9m9' - control mouse
%   - 'b9m8' - noise exposed mouse (96.5 dB)

% % example of 97 dB noise exposed mouse
% MOUSE_LABEL = 'b9m8'; % use lowercase
% MOUSE_DATES = {'20210916', '20210918', '20210925'}; % baseline, post-24 hour, and post-1week time points

% MOUSE_LABEL = 'b9m6'; % use lowercase
% MOUSE_DATES = {'20210913', '20210918', '20210926'}; % baseline, post-24 hour, and post-1week time points

% % % example of 99 dB noise exposed mouse
MOUSE_LABEL = 'b8m6'; % use lowercase
MOUSE_DATES = {'20210911', '20210917', '20210922'}; % baseline, post-24 hour, and post-1week time points

% example of control mouse littermate of 97 dB noise exposed mouse above
% MOUSE_LABEL = 'b9m9'; % use lowercase
% MOUSE_DATES = {'20210916', '20210920', '20210926'}; % baseline, post-24 hour, and post-1week time points

% MOUSE_LABEL = 'b9m10'; % use lowercase
% MOUSE_DATES = {'20210916', '20210920', '20210925'}; % baseline, post-24 hour, and post-1week time points
% 
% MOUSE_LABEL = 'b4m4'; % use lowercase
% MOUSE_DATES = {'20210722', '20210728', '20210803'}; % baseline, post-24 hour, and post-1week time points

% INITIALIZE CONSTANTS FOR PLOTTING ABR ANALYSIS
BATCHGROUPS = {5:6, [2:4,7,9], 8, 2:9}; % group by noise level: 96 dB, 97 dB, and 99 dB respectively.
GROUP_NAMES = {'96 dB', '97 dB', '99 dB', 'all'};
N_GROUPS = numel(BATCHGROUPS);
FREQ = [8000, 16000, 32000];
n_freq = numel(FREQ);
COL_METRICVALS_FIRSTCOL = 25; % first column index of metric values columns in finalTable_metadata table
CONTROL_NOISE_LABEL = {'control', 'noise'};
SIDES = ["left", "right"];
% SIDES = ["both"];
n_sides = numel(SIDES);
SIDES_SYMBOLS = ['x', 'o'];
NAME_TIMEPOINTS = ["Pre", "Post 24h", "Post 1w"];
P_CRIT = 0.05;
LINEWIDTH = 2;
FONTSIZE = 14;

NORMALIZE_ABR = 0; % 9-14-22: binary switch to normalize or not normalize wave 1 amp by maximum amplitude


% METRIC_NAMES = {'Liberman', 'Oghalai', 'Innerprod', 'Innerprod_AUC', 'D_rms', 'D_z'};
% n_metrics = numel(METRIC_NAMES);
% YLIM_METRIC = {[-0.4, 1], [0, 6250], [0, 2.5*10^-9], [0.5, 1]};
% YLABEL_METRIC = {'Correlation coefficient (a.u.)', 'Max p-p amplitude (nV)', 'Mean inner product (nv^2)', 'AUC (a.u.)'}; 

METRIC_NAMES = {'Wave1amp'};
n_metrics = numel(METRIC_NAMES);
YLIM_METRIC = {[0, 8500]};
YLABEL_METRIC = {'Wave 1 amplitude (nV)'}; 
XLIM_METRIC = [-30, 100];
YLIM_METRIC_2 = {[-1.5, 0.7], [-3000, 1600]};
YLABEL_METRIC_2 = {'\Delta Wave 1 amplitude (D'')', '\Delta Wave 1 amplitude (nV)'}; 


% METRIC_NAMES = {'Innerprod_window'};
% n_metrics = numel(METRIC_NAMES);
% YLIM_METRIC = {[-1*10^-12, 1*10^-11]};
% YLABEL_METRIC = {'Inner prod window (nV^2)'}; 

% METRIC_NAMES = {'Wave1amp', 'Wave1amp_auc', 'Innerprod_window', 'Innerprod_window_auc'};
% n_metrics = numel(METRIC_NAMES);
% % YLIM_METRIC = {[0, 6250], [-1*10^-14, 1*10^-13]};
% YLIM_METRIC = {[0, 6250], [0.45, 1], [-5000, 5*10^6], [0.45, 1]};
% YLABEL_METRIC = {'Wave 1 amplitude (nV)', 'Wave 1 amp AUC (a.u.)', 'Inner prod window (nV^2)', 'Inner prod window AUC (a.u.)'}; 

METRIC_THRESHOLD = "Innerprod_auc";

%% Load ABR data. This is a finalTable_metadata variable that was created and saved by save_finalTable_metadata.m
if exist('finalTable_metadata','var') == 0
    load(LOAD_PATH)
end

% %% Load and analyze ABR data - calculate metrics and thresholds
% % Initialize cache variables
% close all
% final_table = [];
%     
% % Iterate over each experiment date
% num_dates = length(MOUSE_DATES);
% for jj=1:num_dates
%     disp(['  Working on file ', num2str(jj), ' out of ', num2str(num_dates), ': ', MOUSE_DATES{jj}])
%     
%     % Get path to ABR data for experimental date
%     filename = [MOUSE_DATES{jj}, '_', MOUSE_LABEL, '_abr_left.csv'];
%     this_path = fullfile(ABR_PATH, MOUSE_DATES{jj}, [MOUSE_DATES{jj}, '_', MOUSE_LABEL], 'analyze');
% 
%     % Analyze ABR
%     output_table = plotstack_averageABR(this_path, filename);
% 
%     % 11-4-21 added
%     % Ensure number of stimulus levels (dB) are same in output table
%     % and cached table results to allow vertical concatenation
%     if ~isempty(final_table)
%         stimulus_levels_output = output_table.Metric_values(1,:);
%         stimulus_levels_cache = final_table.Metric_values(1,:);
% 
%         % check if stimulus levels are same
%         num_levels_output = numel(stimulus_levels_output);
%         num_levels_cache = numel(stimulus_levels_cache);
%         if num_levels_output ~= num_levels_cache
%             % Merge stimulus levels of two tables
%             % cache table
%             outstanding_levels_output = setdiff(stimulus_levels_output, stimulus_levels_cache);
%             if ~isempty(outstanding_levels_output)
%                 % Insert missing levels
%                 add_cols = [];
%                 for z=outstanding_levels_output
%                     new_col = NaN(size(final_table.Metric_values, 1), 1);
%                     new_col(1:3) = z;
%                     add_cols = [add_cols, new_col];
%                 end
%                 new_final_table_metricvals = [final_table.Metric_values, add_cols];
%                 [new_stimulus_levels_cache, ind_levels_sorted] = sort(new_final_table_metricvals(1,:), 'descend');
%                 new_final_table_metricvals = new_final_table_metricvals(:, ind_levels_sorted);
%                 final_table.Metric_values = new_final_table_metricvals;
%             end
% 
%             % output table
%             outstanding_levels_cache = setdiff(stimulus_levels_cache, stimulus_levels_output);
%             if ~isempty(outstanding_levels_cache)
%                 % Insert missing levels
%                 add_cols = [];
%                 for z=outstanding_levels_cache
%                     new_col = NaN(size(output_table.Metric_values, 1), 1);
%                     new_col(1:3) = z;
%                     add_cols = [add_cols, new_col];
%                 end
%                 new_output_table_metricvals = [output_table.Metric_values, add_cols];
%                 [new_stimulus_levels_output, ind_levels_sorted] = sort(new_output_table_metricvals(1,:), 'descend');
%                 new_output_table_metricvals = new_output_table_metricvals(:, ind_levels_sorted);
%                 output_table.Metric_values = new_output_table_metricvals;
%             end
%         end
%     end
% 
%     % Cache table results
%     final_table = vertcat(final_table, output_table);
% end
% 
% % % Write combine output results to excel sheet
% % % [~, baseFileNameNoExt, ~] = fileparts(filename);
% % k = strfind(path, '2021');
% % kk = strfind(path, filesep);
% % is_postslash = kk>k(end);
% % if any(is_postslash)
% %     ind_postslash = kk(is_postslash);
% %     baseFileNameNoExt = path(k:ind_postslash(1)-1);
% % else
% %     baseFileNameNoExt = path(k:end);
% % end
% % 
% % if IS_ABR
% %     file_ext = '_abr_output.xlsx';
% % else
% %     file_ext = '_dpoae_output.xlsx';
% % end
% % table_filename = [baseFileNameNoExt, file_ext];
% % writetable(final_table, fullfile(SAVE_PATH, table_filename), 'Sheet', 1, 'Range', 'A1')
% 
% 
% %% Add metadata to table
% finalTable = final_table;
% n_rows = size(finalTable, 1);
% for k=1:n_rows
%     disp(['Working on metadata row ', num2str(k), ' out of ', num2str(n_rows)])
%     date_name = finalTable.Filenames(k);
%     [date, name, studytype, side, metadata] = get_mousefile_metadata(date_name);
%     metadata.('Date') = date;
%     metadata.('Side') = side;
%     metadata.('Studytype') = studytype;
%     if k==1
%         metadata_aggregate = metadata;
%     else
%         metadata_aggregate = vertcat(metadata_aggregate,metadata);
%     end
% end
% 
% finalTable_metadata = [metadata_aggregate, finalTable];

%% Create filters for selecting data by metric type, noise/control, timepoint
COLS_METRICVALS = COL_METRICVALS_FIRSTCOL:size(finalTable_metadata, 2); % column numbers of metric values columns in finalTable_metadata table

% Metric filters
% is_liberman = contains(finalTable_metadata.('Metric'), 'Liberman', 'IgnoreCase',true);
% is_oghalai = contains(finalTable_metadata.('Metric'), 'Oghalai', 'IgnoreCase',true);
% is_innerprod_auc = contains(finalTable_metadata.('Metric'), 'Innerprod_auc', 'IgnoreCase',true);
% is_innerprod = strcmpi(finalTable_metadata.('Metric'), 'Innerprod');
% metric_selector = {is_liberman, is_oghalai, is_innerprod, is_innerprod_auc};
metric_selector = cell(1, n_metrics);
for i = 1:n_metrics
    metric_selector{i} = strcmpi(finalTable_metadata.('Metric'), METRIC_NAMES{i});
end

% Noise control filter
is_noise = finalTable_metadata.('Is_noise_exposed');
is_control = ~is_noise;
control_noise_selector = {is_control, is_noise};

% Timepoint
is_pre = strcmp(string(finalTable_metadata.('Date_1')), string(finalTable_metadata.('Date')));
is_post24 = strcmp(string(finalTable_metadata.('Date_2')), string(finalTable_metadata.('Date')));
is_post1w = strcmp(string(finalTable_metadata.('Date_3')), string(finalTable_metadata.('Date')));
timepoints = {is_pre, is_post24, is_post1w};
n_timepoints = numel(timepoints);

% 9-13-22: select mouse
is_mouse = strcmpi(finalTable_metadata.Mouse_name, MOUSE_LABEL);

% 10-7-22: threshold metric
threshold_metric_selector = strcmpi(finalTable_metadata.('Metric'), METRIC_THRESHOLD); 

%% Plot level functions for wave 1 amplitude and time window inner product metrics
colors = ['b', 'g', 'r'];

for mm = 1:n_metrics % select metric type
    fig = figure;
    t = tiledlayout(n_sides, n_freq, 'TileIndexing', 'rowmajor');
    
    for nn = 1:n_sides
        if strcmpi(SIDES{nn}, "both")
            is_side = ones(size(finalTable_metadata.Side));
        else
            is_side = finalTable_metadata.Side == SIDES{nn};
        end
        
        n_pts = zeros(n_freq, n_timepoints);
        for fff = 1:n_freq
            % one figure per frequency
            isfreq = finalTable_metadata.('Frequency')==FREQ(fff);

            nexttile(t)
            hold on
            for ii = 1:n_timepoints
                combined_filter = is_mouse & is_side & isfreq & metric_selector{mm} & timepoints{ii};
                data_selected = finalTable_metadata(combined_filter, :);
    %             x = data_selected.Properties.VariableNames(COLS_METRICVALS);
                x = finalTable_metadata(1, COLS_METRICVALS);
                y = data_selected(:, COLS_METRICVALS);

                n_pts(fff, ii) = sum(isfreq);
                ymean = nanmean(y{:,:}, 1);
%                 yste = nanstd(y{:,:}, 1)/sqrt(n_pts(fff, ii));

                % sort points in order of ascending x value to avoid criss
                % crossing of plot lines.
    %             x_double = cellfun(@str2double, x); % convert cell array of char vectors to double matrix
                x_double = x.(1); % convert table 1x1 to double matrix
                [sortedX, sortIndex] = sort(x_double);
                sortedYmean = ymean(sortIndex);

                % Remove nan values to plot uninterrupted curve
                idx = ~(isnan(sortedYmean));
                sortedX_real = sortedX(idx);
                sortedYmean_real = sortedYmean(idx);

                % For wave 1 amp level function, place asterisk over point if
                % wave 1 amplitude distribution of single traces is
                % statistically different from noise distribution. 
                if strcmpi(METRIC_NAMES{mm}, 'Wave1amp')
                    data_selected_dist = finalTable_metadata(combined_filter, :);
                    distribution = data_selected_dist.('Distribution'){1};
                    num_levels = numel(distribution);
                    for aa = 1:num_levels % ascending levels
                        ind = num_levels - aa + 1;
                        is_different_noise = ttest2(distribution{ind}, distribution{end}, 'Tail', 'both');
                        if is_different_noise
                            text(sortedX_real(aa), sortedYmean_real(aa), '*', 'Color', colors(ii), 'FontSize', FONTSIZE, 'FontWeight', 'bold') % asterisk on wave 1 amp curve if significantly greater than noise level's
                        end
                    end

                    % Place red mark above curves if 1 week postexposure level is
                    % significantly distinct from pre exposure level.
                    if ii==n_timepoints
                        data_selected_dist_pre = finalTable_metadata(is_mouse & is_side & isfreq & metric_selector{mm} & is_pre, :);
                        distribution_pre = data_selected_dist_pre.('Distribution'){1};
                        for aa = 1:num_levels % ascending levels
                            ind = num_levels - aa + 1;
                            is_different_pre_post1w = ttest2(distribution{ind}, distribution_pre{ind}, 'Tail', 'both');
                            if is_different_pre_post1w
                                text(sortedX_real(aa), sortedYmean_real(aa)+1500, '*', 'Color', colors(ii), 'FontSize', FONTSIZE*2, 'FontWeight', 'bold') % asterisk above wave 1 amp curve if significantly greater at 1 week c/w baseline
                            end
                        end
                    end
                    
                    % Calculate standard deviation of single trace wave 1
                    % amplitude values to plot their distribution
                    ystd = cellfun(@nanstd, distribution); % Descending order of stimulus levels
%                     sortedYstd = ystd(sortIndex);
                    sortedYstd_real = flip(ystd); % Ascending order of stimulus levels
                end

                % 9-6-22: try normalizing wave1amp level function of single
                % mouse by max stimulus amplitude to see if control mouse then
                % has no change in wave 1 amp at 1 week
                if NORMALIZE_ABR
                    sortedYmean_real = sortedYmean_real/sortedYmean_real(end);
                    
                    % Propagation of errors to adjust standard deviations
                    % of non-max amplitude distributions
                    sortedYstd_real = sqrt(sortedYstd_real.^2 + (sortedYstd_real(end)^2*(sortedYmean_real/sortedYmean_real(end)).^2))/sortedYmean_real(end);
                end
                
                % Plot
                errorbar(sortedX_real, sortedYmean_real, sortedYstd_real, colors(ii), 'LineWidth', LINEWIDTH)
%                 plot(sortedX_real, sortedYmean_real, colors(ii), 'LineWidth', LINEWIDTH)

                
                % For inner product level function, place asterisk above curves if 1 week postexposure level is
                % significantly distinct from pre exposure level.
                if strcmpi(METRIC_NAMES{mm}, 'Innerprod_window') && ii==n_timepoints
                    data_selected_dist = finalTable_metadata(isfreq & metric_selector{mm} & timepoints{ii}, :);
                    distribution = data_selected_dist.('Distribution'){1};
                    data_selected_dist_pre = finalTable_metadata(isfreq & metric_selector{mm} & is_pre, :);
                    distribution_pre = data_selected_dist_pre.('Distribution'){1};
                    num_levels = numel(distribution);
                    for aa = 1:num_levels % ascending levels
                        ind = num_levels - aa + 1;
                        is_different_pre_post1w = ttest2(distribution{ind}, distribution_pre{ind}, 'Tail', 'both');
                        if is_different_pre_post1w
                            text(sortedX_real(aa), sortedYmean_real(aa)+6*10^5, '*', 'Color', colors(ii), 'FontSize', FONTSIZE*2, 'FontWeight', 'bold') % asterisk above inner product curve if significantly greater than noise level's
                        end
                    end
                elseif (strcmpi(METRIC_NAMES{mm}, 'Wave1amp_auc') || strcmpi(METRIC_NAMES{mm}, 'Innerprod_window_auc')) && ii==n_timepoints
                    % If AUC level function, calculate DeLong covariance matrix
                    % to determine 95% confidence intervals and if AUC curves
                    % pre and 1w post exposure are significantly different.

                    % Distributions of innerprods with signal
                    data_selected_dist = finalTable_metadata(isfreq & metric_selector{mm} & timepoints{ii}, :);
                    distribution = data_selected_dist.('Distribution'){1};
                    data_selected_dist_pre = finalTable_metadata(isfreq & metric_selector{mm} & is_pre, :);
                    distribution_pre = data_selected_dist_pre.('Distribution'){1};
                    num_levels = numel(distribution);
                    for aa = 1:num_levels % ascending levels
                        ind = num_levels - aa + 1;

                        signal_pre = distribution_pre{ind}; % signal distribution before noise exposure
                        noise_pre = distribution_pre{end}; % noise distribution before noise exposure
                        signal_post = distribution{ind};
                        noise_post = distribution{end};

                        % Reshape vectors as column vectors
                        signal_pre = reshape(signal_pre, [numel(signal_pre), 1]);
                        noise_pre = reshape(noise_pre, [numel(noise_pre), 1]);
                        signal_post = reshape(signal_post, [numel(signal_post), 1]);
                        noise_post = reshape(noise_post, [numel(noise_post), 1]);

                        % Calculate probability that AUCs of pre and post noise
                        % inner product distrubtions are different (compared
                        % with their respective no signal distributions)
                        data_pre = [signal_pre; noise_pre];
                        data_post1w = [signal_post; noise_post];
                        bm = [data_pre, data_post1w];
                        Disease = [ones(length(signal_pre), 1); zeros(length(noise_pre), 1)];
                        AUC_outputs = AUC_compare_correlated(bm,  Disease); % cell of outputs: {title,AUC1string,AUC2string,AUCdiffstring,prob}
                        prob = AUC_outputs{end};

                        % Plot asterisk above plot point if different
                        is_different_pre_post1w = prob <= P_CRIT ;
                        if is_different_pre_post1w
                            text(sortedX_real(aa), sortedYmean_real(aa), '*', 'Color', colors(ii), 'FontSize', FONTSIZE*2, 'FontWeight', 'bold') % asterisk on plot
                        end
                    end
                end
            end
            hold off
            legend(NAME_TIMEPOINTS, 'Location','northwest')

            xlabel('Amplitude (dB)')
            ylabel(YLABEL_METRIC{mm})
            ylim(YLIM_METRIC{mm})
            if NORMALIZE_ABR
                ylim([-0.1, 1.1]) % For normalized wave1amp
            end
            xlim(XLIM_METRIC)
    %         title([METRIC_NAMES{mm}, ' ', ', ', num2str(FREQ(fff)), ' Hz. N=', num2str(n_pts(fff, 1))]); % TODO: fix n_pts, which is currently 33 times (n_timepts*num_levels) the actual number of mice
            num_single_traces = numel(distribution{1});
            title([num2str(FREQ(fff)), ' Hz, ', SIDES{nn}, ', ', num2str(num_single_traces), ' traces'], 'Interpreter','none');
            set(gca,'FontSize', FONTSIZE) % Set text size of current axis

%             % Save figure to results folder
%             save_filename = ['levelfunction_', METRIC_NAMES{mm}, '_freq', num2str(FREQ(fff)), '_singlemouse_', MOUSE_LABEL, '_', SIDES{nn}];
%             savefig_multiformat(gcf, SAVE_PATH, save_filename)
        end
    end
    
    % Label the tiled layout figure
%     t.TileSpacing = 'none';
    t.TileSpacing = 'tight';
    t.Padding = 'tight';
    if NORMALIZE_ABR
        plot_title = [METRIC_NAMES{mm}, '-', MOUSE_LABEL, '-normalized'];
    else
        plot_title = [METRIC_NAMES{mm}, '-', MOUSE_LABEL];
    end
    title(t, plot_title, 'FontSize', 24)
%         ylabel(t, 'Amplitude (nV)', 'FontSize', 24)
%         xlabel(t, 'Time (ms)', 'FontSize', 24)

    % Maximize figure window size
    fig.WindowState = 'maximized';

    % Save figure to results folder
    save_filename = ['levelfunction_singlemouse_', plot_title];
    savefig_multiformat(gcf, SAVE_PATH, save_filename)
end

%% 9-20-2022 Plot change in wave 1 amplitude from baseline in individual mice, and D' of change
colors = ['b', 'g', 'r'];

for mm = 1:n_metrics % select metric type
    
    for nn = 1:n_sides
        if strcmpi(SIDES{nn}, "both")
            is_side = ones(size(finalTable_metadata.Side));
        else
            is_side = finalTable_metadata.Side == SIDES{nn};
        end
        
        fig = figure;
        t = tiledlayout(2, n_freq, 'TileIndexing', 'columnmajor');
        
        n_pts = zeros(n_freq, n_timepoints);
        for fff = 1:n_freq
            % one figure per frequency
            isfreq = finalTable_metadata.('Frequency')==FREQ(fff);

            for i=1:2 % Plot D' then wave 1 amp change
                nexttile(t)
                hold on
                for ii = 1:n_timepoints
                    combined_filter = is_mouse & is_side & isfreq & metric_selector{mm} & timepoints{ii};
                    data_selected = finalTable_metadata(combined_filter, :);
        %             x = data_selected.Properties.VariableNames(COLS_METRICVALS);
                    x = finalTable_metadata(1, COLS_METRICVALS);
                    y = data_selected(:, COLS_METRICVALS);

                    n_pts(fff, ii) = sum(isfreq);
                    ymean = nanmean(y{:,:}, 1);
    %                 yste = nanstd(y{:,:}, 1)/sqrt(n_pts(fff, ii));

                    % sort points in order of ascending x value to avoid criss
                    % crossing of plot lines.
        %             x_double = cellfun(@str2double, x); % convert cell array of char vectors to double matrix
                    x_double = x.(1); % convert table 1x1 to double matrix
                    [sortedX, sortIndex] = sort(x_double);
                    sortedYmean = ymean(sortIndex);

                    % Remove nan values to plot uninterrupted curve
                    idx = ~(isnan(sortedYmean));
                    sortedX_real = sortedX(idx);
                    sortedYmean_real = sortedYmean(idx);

                    % Calculate change in wave 1 amplitude from baseline to
                    % post-noise date
                    d_prime = zeros(num_levels, 1);
                    delta_wave1amp = zeros(num_levels, 1);
                    ste_delta_wave1amp = zeros(num_levels, 1);
                    if strcmpi(METRIC_NAMES{mm}, 'Wave1amp')
                        num_levels = numel(distribution);

                        if ii==1
                            distribution_pre = data_selected.('Distribution'){1};
                        else                        
                            for aa = 1:num_levels % ascending levels
                                ind = num_levels - aa + 1;
                                distribution_post = data_selected.('Distribution'){1};

                                this_pre_wave1amp_dist = distribution_pre{ind};
                                this_post_wave1amp_dist = distribution_post{ind};

                                % Calculate D' for this stimulus level
                                n1 = numel(this_pre_wave1amp_dist);
                                mean1 = mean(this_pre_wave1amp_dist);
                                std1 = std(this_pre_wave1amp_dist);
                                n2 = numel(this_post_wave1amp_dist);
                                mean2 = mean(this_post_wave1amp_dist);
                                std2 = std(this_post_wave1amp_dist);
                                [npool,meanpool,stdpool] = pooledmeanstd(n1,mean1,std1,n2,mean2,std2);
                                d_prime(aa) = (mean2 - mean1)/stdpool;

                                delta_wave1amp(aa) = mean2 - mean1;
                                ste1 = std1/sqrt(n1);
                                ste2 = std2/sqrt(n2);
                                ste_delta_wave1amp(aa) = sqrt(ste1^2 + ste2^2);
    %                             is_different_pre_post1w = ttest2(distribution{ind}, distribution_pre{ind}, 'Tail', 'both');
    %                             if is_different_pre_post1w
    %                                 text(sortedX_real(aa), sortedYmean_real(aa)+1500, '*', 'Color', colors(ii), 'FontSize', FONTSIZE*2, 'FontWeight', 'bold') % asterisk above wave 1 amp curve if significantly greater at 1 week c/w baseline
    %                             end
                            end
                        end

                        % Calculate standard deviation of single trace wave 1
                        % amplitude values to plot their distribution
                        ystd = cellfun(@nanstd, distribution); % Descending order of stimulus levels
    %                     sortedYstd = ystd(sortIndex);
                        sortedYstd_real = flip(ystd); % Ascending order of stimulus levels
                    end

                    % Remove stimulus levels that are not divisible by 10
                    isdiv10 = mod(sortedX_real, 10)==0;
                    sortedX_real_isdiv10 = sortedX_real(isdiv10);
                    delta_wave1amp_isdiv10 = delta_wave1amp(isdiv10);
                    ste_delta_wave1amp_isdiv10 = ste_delta_wave1amp(isdiv10);
                    
                    if i==1 % Plot D'
                        plot(sortedX_real(isdiv10), d_prime(isdiv10), colors(ii), 'LineWidth', LINEWIDTH)
                    else
                        
                        % Plot change in wave 1 amplitude from baseline to after
                        % noise exposure
                        errorbar(sortedX_real_isdiv10, delta_wave1amp_isdiv10, ste_delta_wave1amp_isdiv10, colors(ii), 'LineWidth', LINEWIDTH);
                        
                        % 10-7-22: Plot vertical line to show threshold
                        combined_filter_threshold = is_mouse & is_side & isfreq & threshold_metric_selector & timepoints{ii};
                        threshold_selected = finalTable_metadata.('Threshold')(combined_filter_threshold);
                        mean_threshold_ears = nanmean(threshold_selected); % in case multiple thresholds are measured for both ears
                        p = xline(mean_threshold_ears, '--', 'Color', colors(ii));
                        set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Do not include vertical threshold lines in legend
                        
                        % Check statistical significance from zero
                        z_vector = delta_wave1amp_isdiv10./ste_delta_wave1amp_isdiv10;
                        p_twoTailed = 2*(1 - normcdf(abs(z_vector)));
                        is_significant = p_twoTailed < P_CRIT;
                        zz = find(is_significant);
                        text(sortedX_real_isdiv10(zz), delta_wave1amp_isdiv10(zz), '*', ...
                            'Color', colors(ii), 'FontSize', FONTSIZE*2, 'FontWeight', 'bold') 
                    end

                end
                
                hold off
                legend(NAME_TIMEPOINTS, 'Location', 'southwest')

                xlabel('Amplitude (dB)')
                ylabel(YLABEL_METRIC_2{i})
                ylim(YLIM_METRIC_2{i})
                xlim(XLIM_METRIC)
                num_single_traces = numel(this_post_wave1amp_dist);
                if i==1
                    title([num2str(FREQ(fff)), ' Hz, D'''], 'Interpreter','none');
                else
                    title([num2str(FREQ(fff)), ' Hz, change in wave 1 amp'], 'Interpreter','none');
                end
                set(gca,'FontSize', FONTSIZE) % Set text size of current axis

    %             % Save figure to results folder
    %             save_filename = ['levelfunction_', METRIC_NAMES{mm}, '_freq', num2str(FREQ(fff)), '_singlemouse_', MOUSE_LABEL, '_', SIDES{nn}];
    %             savefig_multiformat(gcf, SAVE_PATH, save_filename)
            end
        end
        
         % Label the tiled layout figure
    %     t.TileSpacing = 'none';
        t.TileSpacing = 'tight';
        t.Padding = 'tight';
        plot_title = [METRIC_NAMES{mm}, ', ', SIDES{nn}, ', ', MOUSE_LABEL, ', ', num2str(num_single_traces), ' traces'];
        title(t, plot_title, 'FontSize', 24)
    %         ylabel(t, 'Amplitude (nV)', 'FontSize', 24)
    %         xlabel(t, 'Time (ms)', 'FontSize', 24)

        % Maximize figure window size
        fig.WindowState = 'maximized';

        % Save figure to results folder
        save_filename = ['levelfunction_singlemouse_change_', plot_title];
        savefig_multiformat(gcf, SAVE_PATH, save_filename)
    end
    
   
end

disp('Done.')