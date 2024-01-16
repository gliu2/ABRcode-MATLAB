% Plot effect size D' for change in group average wave 1 amplitude
% post24hour and post1week compared with pre (baseline).
%
% Group mice by noise exposure group
% Code adapted from plot_ABRthreshold_levelfunctions_onemouse.m
%
% SAVE_PATH automatically set to "Results" folder to save output files
%
% Dependencies: plotstack_averageABR.m, pooledmeanstd.m, ttest2_mean_ste
%
% 9/14/2022 George Liu

close all
opengl('save', 'software') % prevent crashing due to low-level graphics error using Sarah Office computer's graphics card

%% Constants
IS_ABR = 1; % MANUALLY SET TO 1 (ABR) OR 0 (DPOAE)
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures
LOAD_PATH = 'd:\users\admin\Documents\George\finalTable_metadata_4-18-22\finalTable_metadata_9-13-22.mat'; % load variable finalTable_metadata

% INITIALIZE CONSTANTS FOR PLOTTING ABR ANALYSIS
BATCHGROUPS = {5:6, [2:4,7,9], 8, 2:9}; % group by noise level: 96 dB, 97 dB, 99 dB, and all respectively.
GROUP_NAMES = {'96 dB', '97 dB', '99 dB', 'all'};
N_GROUPS = numel(BATCHGROUPS);
FREQ = [8000, 16000, 32000];
n_freq = numel(FREQ);
COL_METRICVALS_FIRSTCOL = 25; % first column index of metric values columns in finalTable_metadata table
CONTROL_NOISE_LABEL = {'control', 'noise'};
% SIDES = ["left", "right"];
SIDES = ["both"];
n_sides = numel(SIDES);
SIDES_SYMBOLS = ['x', 'o'];
NAME_TIMEPOINTS = ["Pre", "Post 24h", "Post 1w"];
P_CRIT = 0.05;
LINEWIDTH = 2;
FONTSIZE = 14;

METRIC_NAMES = {'Wave1amp'};
n_metrics = numel(METRIC_NAMES);
YLIM_METRIC = {[-2.1, 2.1]};
YLABEL_METRIC = {'\Delta Wave 1 amplitude (D'')'}; 
YLIM_METRIC_RAW = {[-3000, 2000]};
YLABEL_METRIC_RAWCHANGE = {'\Delta Wave 1 amplitude (nV)'}; 

%% Load ABR data. This is a finalTable_metadata variable that was created and saved by save_finalTable_metadata.m
if exist('finalTable_metadata','var') == 0
    load(LOAD_PATH)
end

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
n_conditions = numel(control_noise_selector);

% Timepoint
is_pre = strcmp(string(finalTable_metadata.('Date_1')), string(finalTable_metadata.('Date')));
is_post24 = strcmp(string(finalTable_metadata.('Date_2')), string(finalTable_metadata.('Date')));
is_post1w = strcmp(string(finalTable_metadata.('Date_3')), string(finalTable_metadata.('Date')));
timepoints = {is_pre, is_post24, is_post1w};
n_timepoints = numel(timepoints);

%% Plot level functions for wave 1 amplitude and time window inner product metrics
colors = ['b', 'g', 'r'];

for mm = 1:n_metrics % Iterate over each metric 
    for gg = 1:N_GROUPS % Iterate over noise group
        is_ingroup = any(finalTable_metadata.Batch == BATCHGROUPS{gg}, 2); % Group filter
        for nn = 1:n_sides
            if strcmpi(SIDES{nn}, "both")
                is_side = ones(size(finalTable_metadata.Side));
            else
                is_side = finalTable_metadata.Side == SIDES{nn};
            end
            
            fig = figure;
            t = tiledlayout(2, n_freq, 'TileIndexing', 'columnmajor');
            
            for fff = 1:n_freq
                % plot each frequency separately
                isfreq = finalTable_metadata.('Frequency')==FREQ(fff);
            
                % Plot effect size (D') of wave 1 amplitude change
                nexttile(t)
                n_pts = zeros(n_conditions, 1);
                handles = [];
                handle_count = 1;
                for cc = 1:n_conditions % Iterate over control then noise mice
                    
                    hold on
                    y_pre = [];
                    for ii = 1:n_timepoints
                        combined_filter = is_side & is_ingroup & control_noise_selector{cc} & isfreq & metric_selector{mm} & timepoints{ii};
                        data_selected = finalTable_metadata(combined_filter, :);
            %             x = data_selected.Properties.VariableNames(COLS_METRICVALS);
                        x = finalTable_metadata(1, COLS_METRICVALS);
                        y = data_selected(:, COLS_METRICVALS);

                        n_pts(fff, ii) = sum(combined_filter);
                        ymean = nanmean(y{:,:}, 1);
        %                 yste = nanstd(y{:,:}, 1)/sqrt(n_pts(gg, fff, ii)); % standard error
                        yste = nanstd(y{:,:}, 1); % Standard deviation



                        % sort points in order of ascending x value to avoid criss
                        % crossing of plot lines.
            %             x_double = cellfun(@str2double, x); % convert cell array of char vectors to double matrix
                        x_double = x.(1); % convert table 1x1 to double matrix
                        [sortedX, sortIndex] = sort(x_double);
                        sortedYmean = ymean(sortIndex);
                        sortedYstd = yste(sortIndex);

                        % Remove nan values to plot uninterrupted curve
                        idx = ~(isnan(sortedYmean));
                        sortedX_real = sortedX(idx);
                        sortedYmean_real = sortedYmean(idx);
                        sortedYstd_real = sortedYstd(idx);
                        num_levels = numel(sortedX);
                       
%                         errorbar(sortedX_real(isdiv10), sortedYmean_real(isdiv10), sortedYstd_real(isdiv10), colors(ii), 'LineWidth', LINEWIDTH)
        %                 plot(sortedX_real, sortedYmean_real, colors(ii), 'LineWidth', LINEWIDTH)

                        % Determine if pre and post are significantly different
                        % at each stimulus level
                        if ii==1
                            y_pre = y{:,:};
                            y_pre = y_pre(:, sortIndex);
                        else
                            y_post = y{:,:};
                            y_post = y_post(:, sortIndex);
                            d_prime = zeros(num_levels, 1);
                            for zz = 1:num_levels % iterate thru stimulus levels
                                this_prey = y_pre(:, zz);
                                this_posty = y_post(:, zz);

                                if all(isnan(this_prey)) || all(isnan(this_posty))
                                    continue
                                end

                                % Remove NaN values
                                this_prey = this_prey(~isnan(this_prey));
                                this_posty = this_posty(~isnan(this_posty));
                                
                                % Calculate D' for this stimulus level
                                n1 = numel(this_prey);
                                mean1 = mean(this_prey);
                                std1 = std(this_prey);
                                n2 = numel(this_posty);
                                mean2 = mean(this_posty);
                                std2 = std(this_posty);
                                [npool,meanpool,stdpool] = pooledmeanstd(n1,mean1,std1,n2,mean2,std2);
                                d_prime(zz) = (mean2 - mean1)/stdpool;
                                
%                                 % 9-16-22 Also calculate change in wave 1
%                                 % amplitude compared to baseline
%                                 delta_wave1amp = mean2 - mean1;
%                                 ste1 = std1/sqrt(n1);
%                                 ste2 = std2/sqrt(n2);
%                                 ste_delta_wave1amp = sqrt(ste1^2 + ste2^2);

                                % Compare with t-test
                                is_different = ttest2(this_prey, this_posty, 'Tail', 'both');
                                if is_different && mod(sortedX(zz), 10)==0
                                    % asterisk on wave 1 amp curve if significantly different from baseline
                                    text(sortedX(zz), d_prime(zz), '*', ...
                                        'Color', colors(ii), 'FontSize', FONTSIZE*2, 'FontWeight', 'bold') 
                                end
                            end
                            
                            % Remove stimulus levels that are not divisible by 10
                            isdiv10 = mod(sortedX, 10)==0;
                            
                            % Plot D' for this post-exposure date
                            if strcmpi(CONTROL_NOISE_LABEL{cc}, 'control') % control lines are dashed
                                handles(handle_count) = plot(sortedX(isdiv10), d_prime(isdiv10), colors(ii), ...
                                    'LineWidth', LINEWIDTH, 'LineStyle', "--");
                            else % noise exposed mouse lines are solid
                                handles(handle_count) = plot(sortedX(isdiv10), d_prime(isdiv10), colors(ii), ...
                                    'LineWidth', LINEWIDTH);
                            end
                            
                            handle_count = handle_count + 1;
                        end

                    end

                    %TODO: add statistical significance test between pre and post
                    %noise curves

                    hold off
                    legend(handles, [strcat(NAME_TIMEPOINTS(2:end), '-control'), strcat(NAME_TIMEPOINTS(2:end), '-noise')], 'Location','northwest', 'AutoUpdate','off')
                    yline(0, ':')
                    yline(-1, ':')
                    yline(1, ':')

                    xlabel('Amplitude (dB)')
                    ylabel(YLABEL_METRIC{mm})
                    ylim(YLIM_METRIC{mm})
        %             ylim([-0.1, 1.1]) % For normalized wave1amp
                    xlim([-30, 100])
            %         title([METRIC_NAMES{mm}, ' ', ', ', num2str(FREQ(fff)), ' Hz. N=', num2str(n_pts(fff, 1))]); % TODO: fix n_pts, which is currently 33 times (n_timepts*num_levels) the actual number of mice
                    title([num2str(FREQ(fff)), ' Hz. N=', num2str(n_pts(fff, ii)), ' ears'], 'Interpreter','none');
                    set(gca,'FontSize', FONTSIZE) % Set text size of current axis
                end
                
                % 9-16-22: Plot raw value of wave 1 amplitude change for group
                % average
                nexttile(t)
                n_pts = zeros(n_conditions, 1);
                handles = [];
                handle_count = 1;
                
                % Initialize control condition wave 1 amplitude change
                % statistics for comparison with noise condition
                control_n = zeros(n_timepoints - 1, num_levels);
                control_mean = zeros(n_timepoints - 1, num_levels);
                control_ste = zeros(n_timepoints - 1, num_levels);
                for cc = 1:n_conditions % Iterate over control then noise mice
                    
                    hold on
                    y_pre = [];
                    for ii = 1:n_timepoints
                        combined_filter = is_side & is_ingroup & control_noise_selector{cc} & isfreq & metric_selector{mm} & timepoints{ii};
                        data_selected = finalTable_metadata(combined_filter, :);
            %             x = data_selected.Properties.VariableNames(COLS_METRICVALS);
                        x = finalTable_metadata(1, COLS_METRICVALS);
                        y = data_selected(:, COLS_METRICVALS);

                        n_pts(fff, ii) = sum(combined_filter);
                        ymean = nanmean(y{:,:}, 1);
        %                 yste = nanstd(y{:,:}, 1)/sqrt(n_pts(gg, fff, ii)); % standard error
                        yste = nanstd(y{:,:}, 1); % Standard deviation



                        % sort points in order of ascending x value to avoid criss
                        % crossing of plot lines.
            %             x_double = cellfun(@str2double, x); % convert cell array of char vectors to double matrix
                        x_double = x.(1); % convert table 1x1 to double matrix
                        [sortedX, sortIndex] = sort(x_double);
                        sortedYmean = ymean(sortIndex);
                        sortedYstd = yste(sortIndex);

                        % Remove nan values to plot uninterrupted curve
                        idx = ~(isnan(sortedYmean));
                        sortedX_real = sortedX(idx);
                        sortedYmean_real = sortedYmean(idx);
                        sortedYstd_real = sortedYstd(idx);
                        num_levels = numel(sortedX);
                       
%                         errorbar(sortedX_real(isdiv10), sortedYmean_real(isdiv10), sortedYstd_real(isdiv10), colors(ii), 'LineWidth', LINEWIDTH)
        %                 plot(sortedX_real, sortedYmean_real, colors(ii), 'LineWidth', LINEWIDTH)

                        % Determine if pre and post are significantly different
                        % at each stimulus level
                        if ii==1
                            y_pre = y{:,:};
                            y_pre = y_pre(:, sortIndex);
                        else
                            y_post = y{:,:};
                            y_post = y_post(:, sortIndex);
                            d_prime = zeros(num_levels, 1);
                            delta_wave1amp = zeros(num_levels, 1);
                            ste_delta_wave1amp = zeros(num_levels, 1);
                            for zz = 1:num_levels % iterate thru stimulus levels
                                this_prey = y_pre(:, zz);
                                this_posty = y_post(:, zz);

                                if all(isnan(this_prey)) || all(isnan(this_posty))
                                    continue
                                end

                                % Remove NaN values
                                this_prey = this_prey(~isnan(this_prey));
                                this_posty = this_posty(~isnan(this_posty));
                                
                                % Calculate D' for this stimulus level
                                n1 = numel(this_prey);
                                mean1 = mean(this_prey);
                                std1 = std(this_prey);
                                n2 = numel(this_posty);
                                mean2 = mean(this_posty);
                                std2 = std(this_posty);
                                [npool,meanpool,stdpool] = pooledmeanstd(n1,mean1,std1,n2,mean2,std2);
                                d_prime(zz) = (mean2 - mean1)/stdpool;
                                
                                % 9-16-22 Also calculate change in wave 1
                                % amplitude compared to baseline
                                delta_wave1amp(zz) = mean2 - mean1;
                                ste1 = std1/sqrt(n1);
                                ste2 = std2/sqrt(n2);
                                ste_delta_wave1amp(zz) = sqrt(ste1^2 + ste2^2);

                                % Compare wave 1 amplitude change between control and noise conditions using 2-sample unpaired t-test
                                if strcmpi(CONTROL_NOISE_LABEL{cc}, 'control') % control condition - store wave 1 amplitude changes
                                    control_n(ii - 1, zz) = n2;
                                    control_mean(ii - 1, zz) = delta_wave1amp(zz);
                                    control_ste(ii - 1, zz) = ste_delta_wave1amp(zz);
                                else  % noise condition - perform unpaired 2 sample t-test to determine if different from control condition
                                    is_different = ttest2_mean_ste(control_n(ii - 1, zz), control_mean(ii - 1, zz),...
                                        control_ste(ii - 1, zz), n2, delta_wave1amp(zz), ste_delta_wave1amp(zz));
                                    if is_different < P_CRIT && mod(sortedX(zz), 10)==0
                                        % asterisk on noise condition's
                                        % delta wave 1 amp curve if
                                        % significantly different from
                                        % control change
                                        text(sortedX(zz), delta_wave1amp(zz), '*', ...
                                            'Color', colors(ii), 'FontSize', FONTSIZE*2, 'FontWeight', 'bold') 
                                    end
                                end
                            end
                            
                            % Remove stimulus levels that are not divisible by 10
                            isdiv10 = mod(sortedX, 10)==0;
                            
                            % Plot raw wave 1 amplitude change for this post-exposure date
                            if strcmpi(CONTROL_NOISE_LABEL{cc}, 'control') % control lines are dashed
                                handles(handle_count) = errorbar(sortedX(isdiv10), delta_wave1amp(isdiv10), ...
                                    ste_delta_wave1amp(isdiv10), colors(ii), 'LineWidth', LINEWIDTH, 'LineStyle', "--");
                            else % noise exposed mouse lines are solid
                                handles(handle_count) = errorbar(sortedX(isdiv10), delta_wave1amp(isdiv10), ...
                                    ste_delta_wave1amp(isdiv10), colors(ii), 'LineWidth', LINEWIDTH);
                            end
                            
                            handle_count = handle_count + 1;
                        end

                    end

                    %TODO: add statistical significance test between pre and post
                    %noise curves

                    hold off
                    legend(handles, [strcat(NAME_TIMEPOINTS(2:end), '-control'), strcat(NAME_TIMEPOINTS(2:end), '-noise')], 'Location','northwest', 'AutoUpdate','off')
                    yline(0, ':')
                    yline(-1, ':')
                    yline(1, ':')

                    xlabel('Amplitude (dB)')
                    ylabel(YLABEL_METRIC_RAWCHANGE{mm})
                    ylim(YLIM_METRIC_RAW{mm})
        %             ylim([-0.1, 1.1]) % For normalized wave1amp
                    xlim([-30, 100])
            %         title([METRIC_NAMES{mm}, ' ', ', ', num2str(FREQ(fff)), ' Hz. N=', num2str(n_pts(fff, 1))]); % TODO: fix n_pts, which is currently 33 times (n_timepts*num_levels) the actual number of mice
%                     title([CONTROL_NOISE_LABEL{cc}, ', ', num2str(FREQ(fff)), ' Hz. N=', num2str(n_pts(fff, ii)), ' ears'], 'Interpreter','none');
                    set(gca,'FontSize', FONTSIZE) % Set text size of current axis
                end
            end

            % Label the tiled layout figure
    %         t.TileSpacing = 'none';
            t.TileSpacing = 'tight';
            t.Padding = 'tight';
            plot_title = [GROUP_NAMES{gg}, ', ', SIDES{nn}, ' side, ', METRIC_NAMES{mm}];
            title(t, plot_title, 'FontSize', 24)
    %         ylabel(t, 'Amplitude (nV)', 'FontSize', 24)
    %         xlabel(t, 'Time (ms)', 'FontSize', 24)

            % Maximize figure window size
            fig.WindowState = 'maximized';

            % Save figure to results folder
            save_filename = ['levelfunction_group', GROUP_NAMES{gg}, '_', SIDES{nn}, '_', METRIC_NAMES{mm}];
            savefig_multiformat(gcf, SAVE_PATH, save_filename)
        end
    end
end

% %% 9-21-2022 Plot distribution of individual mouse D' and wave 1 amp change for one frequency and stimulus level
% %% TODO: CODE IS INCOMPLETE IN THIS SECTION 
% colors = ['b', 'g', 'r'];
% THIS_FREQUENCY = 32000;
% THIS_LEVEL = 60;
% 
% isfreq = finalTable_metadata.('Frequency')==THIS_FREQUENCY;
% 
% for mm = 1:n_metrics % Iterate over each metric 
%     for nn = 1:n_sides
%         if strcmpi(SIDES{nn}, "both")
%             is_side = ones(size(finalTable_metadata.Side));
%         else
%             is_side = finalTable_metadata.Side == SIDES{nn};
%         end
%     
%         for gg = 1:N_GROUPS % Iterate over noise group
%             is_ingroup = any(finalTable_metadata.Batch == BATCHGROUPS{gg}, 2); % Group filter
%         
%             fig = figure;
%             t = tiledlayout(4, N_GROUPS, 'TileIndexing', 'columnmajor');
%             
%             dprime_control_noise = cell(n_conditions, 1);
%             delta_wave1amp_control_noise = cell(n_conditions, 1);
%             ste_delta_wave1amp_control_noise = cell(n_conditions, 1);
%             
%             for cc = 1:n_conditions % Iterate over control then noise mice
%                 mouse_filter = is_ingroup & control_noise_selector{cc};
%                 data_selected = finalTable_metadata(mouse_filter, :);
%                 these_mice = unique(data_selected.('Mouse_name'));
%                 num_mice = numel(these_mice);
%             
%                 d_prime_pre = zeros(num_mice, 1);
%                 d_prime_post24h = zeros(num_mice, 1);
%                 d_prime_post1w = zeros(num_mice, 1);
%                 delta_wave1amp_pre = zeros(num_mice, 1);
%                 delta_wave1amp_post24h = zeros(num_mice, 1);
%                 delta_wave1amp_post1w = zeros(num_mice, 1);
%                 ste_delta_wave1amp_pre = zeros(num_mice, 1);
%                 ste_delta_wave1amp_post24h = zeros(num_mice, 1);
%                 ste_delta_wave1amp_post1w = zeros(num_mice, 1);
%                 
%                 for i=1:num_mice
%                     this_mouse = these_mice{i};
%                     is_mouse = strcmpi(finalTable_metadata.Mouse_name, this_mouse);
%             
%                     [d_prime_cache, delta_wave1amp_cache, ste_delta_wave1amp_cache, ...
%                         sortedX_real_cache, sortedX_real_significantchange_cache, ...
%                         delta_wave1amp_significantchange_cache] = get_wave1amp_change_singlemouse(finalTable_metadata, ...
%                         is_mouse, is_side, isfreq, metric_selector{mm}, timepoints);
%                     
%                     IS_THIS_LEVEL = sortedX_real_cache{1}==THIS_LEVEL;
%                     
%                     % Cache results
%                     d_prime_pre(i) = d_prime_cache{1}(IS_THIS_LEVEL);
%                     d_prime_post24h(i) = d_prime_cache{2}(IS_THIS_LEVEL);
%                     d_prime_post1w(i) = d_prime_cache{3}(IS_THIS_LEVEL);
%                     delta_wave1amp_pre(i) = delta_wave1amp_cache{1}(IS_THIS_LEVEL);
%                     delta_wave1amp_post24h(i) = delta_wave1amp_cache{2}(IS_THIS_LEVEL);
%                     delta_wave1amp_post1w(i) = delta_wave1amp_cache{3}(IS_THIS_LEVEL);
%                     ste_delta_wave1amp_pre(i) = ste_delta_wave1amp_cache{1}(IS_THIS_LEVEL);
%                     ste_delta_wave1amp_post24h(i) = ste_delta_wave1amp_cache{2}(IS_THIS_LEVEL);
%                     ste_delta_wave1amp_post1w(i) = ste_delta_wave1amp_cache{3}(IS_THIS_LEVEL);
%                     
%                 end
%                 
%                 dprime_control_noise{cc} = {d_prime_pre, d_prime_post24h, d_prime_post1w};
%                 delta_wave1amp_control_noise{cc} = {delta_wave1amp_pre, delta_wave1amp_post24h, delta_wave1amp_post1w};
%                 ste_delta_wave1amp_control_noise{cc} = {ste_delta_wave1amp_pre, ste_delta_wave1amp_post24h, ste_delta_wave1amp_post1w};
%             end
%             
%             % Plot D' and wave 1 amp change individual mouse distributions
%             % at 24h and 1w post noise exposure, for noise and control mice
%             % for each group
%             
%             % D' at 24h and 1 w
%             box_labels = {'24h', '1w', 'control'};
%             for ww = 2:3
%                 nexttile(t)
%                 x = [dprime_control_noise{2}{ww}, dprime_control_noise{1}{ww}]; % [noise, control]
%                 box_labels = {time_labels{ww}, 'control'};
%                 boxplot(x, 'Labels', box_labels)
%                 ylabel('D''')
%             end
%             
%             % Wave 1 amp change at 24h and 1 w
%             for ww = 2:3
%                 nexttile(t)
%                 x = [delta_wave1amp_control_noise{2}{ww}, delta_wave1amp_control_noise{1}{ww}]; % [noise, control]
%                 box_labels = {time_labels{ww}, 'control'};
%                 boxplot(x, 'Labels', box_labels)
%             end
% 
%             
%                 % Plot effect size (D') of wave 1 amplitude change
%                 nexttile(t)
%                 n_pts = zeros(n_conditions, 1);
%                 handles = [];
%                 handle_count = 1;
%                 
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                     hold on
%                     y_pre = [];
%                     for ii = 1:n_timepoints
%                         combined_filter = is_side & is_ingroup & control_noise_selector{cc} & isfreq & metric_selector{mm} & timepoints{ii};
%                         data_selected = finalTable_metadata(combined_filter, :);
%             %             x = data_selected.Properties.VariableNames(COLS_METRICVALS);
%                         x = finalTable_metadata(1, COLS_METRICVALS);
%                         y = data_selected(:, COLS_METRICVALS);
% 
%                         n_pts(fff, ii) = sum(combined_filter);
%                         ymean = nanmean(y{:,:}, 1);
%         %                 yste = nanstd(y{:,:}, 1)/sqrt(n_pts(gg, fff, ii)); % standard error
%                         yste = nanstd(y{:,:}, 1); % Standard deviation
% 
% 
% 
%                         % sort points in order of ascending x value to avoid criss
%                         % crossing of plot lines.
%             %             x_double = cellfun(@str2double, x); % convert cell array of char vectors to double matrix
%                         x_double = x.(1); % convert table 1x1 to double matrix
%                         [sortedX, sortIndex] = sort(x_double);
%                         sortedYmean = ymean(sortIndex);
%                         sortedYstd = yste(sortIndex);
% 
%                         % Remove nan values to plot uninterrupted curve
%                         idx = ~(isnan(sortedYmean));
%                         sortedX_real = sortedX(idx);
%                         sortedYmean_real = sortedYmean(idx);
%                         sortedYstd_real = sortedYstd(idx);
%                         num_levels = numel(sortedX);
%                        
% %                         errorbar(sortedX_real(isdiv10), sortedYmean_real(isdiv10), sortedYstd_real(isdiv10), colors(ii), 'LineWidth', LINEWIDTH)
%         %                 plot(sortedX_real, sortedYmean_real, colors(ii), 'LineWidth', LINEWIDTH)
% 
%                         % Determine if pre and post are significantly different
%                         % at each stimulus level
%                         if ii==1
%                             y_pre = y{:,:};
%                             y_pre = y_pre(:, sortIndex);
%                         else
%                             y_post = y{:,:};
%                             y_post = y_post(:, sortIndex);
%                             d_prime = zeros(num_levels, 1);
%                             for zz = 1:num_levels % iterate thru stimulus levels
%                                 this_prey = y_pre(:, zz);
%                                 this_posty = y_post(:, zz);
% 
%                                 if all(isnan(this_prey)) || all(isnan(this_posty))
%                                     continue
%                                 end
% 
%                                 % Remove NaN values
%                                 this_prey = this_prey(~isnan(this_prey));
%                                 this_posty = this_posty(~isnan(this_posty));
%                                 
%                                 % Calculate D' for this stimulus level
%                                 n1 = numel(this_prey);
%                                 mean1 = mean(this_prey);
%                                 std1 = std(this_prey);
%                                 n2 = numel(this_posty);
%                                 mean2 = mean(this_posty);
%                                 std2 = std(this_posty);
%                                 [npool,meanpool,stdpool] = pooledmeanstd(n1,mean1,std1,n2,mean2,std2);
%                                 d_prime(zz) = (mean2 - mean1)/stdpool;
%                                 
% %                                 % 9-16-22 Also calculate change in wave 1
% %                                 % amplitude compared to baseline
% %                                 delta_wave1amp = mean2 - mean1;
% %                                 ste1 = std1/sqrt(n1);
% %                                 ste2 = std2/sqrt(n2);
% %                                 ste_delta_wave1amp = sqrt(ste1^2 + ste2^2);
% 
%                                 % Compare with t-test
%                                 is_different = ttest2(this_prey, this_posty, 'Tail', 'both');
%                                 if is_different && mod(sortedX(zz), 10)==0
%                                     % asterisk on wave 1 amp curve if significantly different from baseline
%                                     text(sortedX(zz), d_prime(zz), '*', ...
%                                         'Color', colors(ii), 'FontSize', FONTSIZE*2, 'FontWeight', 'bold') 
%                                 end
%                             end
%                             
%                             % Remove stimulus levels that are not divisible by 10
%                             isdiv10 = mod(sortedX, 10)==0;
%                             
%                             % Plot D' for this post-exposure date
%                             if strcmpi(CONTROL_NOISE_LABEL{cc}, 'control') % control lines are dashed
%                                 handles(handle_count) = plot(sortedX(isdiv10), d_prime(isdiv10), colors(ii), ...
%                                     'LineWidth', LINEWIDTH, 'LineStyle', "--");
%                             else % noise exposed mouse lines are solid
%                                 handles(handle_count) = plot(sortedX(isdiv10), d_prime(isdiv10), colors(ii), ...
%                                     'LineWidth', LINEWIDTH);
%                             end
%                             
%                             handle_count = handle_count + 1;
%                         end
% 
%                     end
% 
%                     %TODO: add statistical significance test between pre and post
%                     %noise curves
% 
%                     hold off
%                     legend(handles, [strcat(NAME_TIMEPOINTS(2:end), '-control'), strcat(NAME_TIMEPOINTS(2:end), '-noise')], 'Location','northwest', 'AutoUpdate','off')
%                     yline(0, ':')
%                     yline(-1, ':')
%                     yline(1, ':')
% 
%                     xlabel('Amplitude (dB)')
%                     ylabel(YLABEL_METRIC{mm})
%                     ylim(YLIM_METRIC{mm})
%         %             ylim([-0.1, 1.1]) % For normalized wave1amp
%                     xlim([-30, 100])
%             %         title([METRIC_NAMES{mm}, ' ', ', ', num2str(FREQ(fff)), ' Hz. N=', num2str(n_pts(fff, 1))]); % TODO: fix n_pts, which is currently 33 times (n_timepts*num_levels) the actual number of mice
%                     title([num2str(FREQ(fff)), ' Hz. N=', num2str(n_pts(fff, ii)), ' ears'], 'Interpreter','none');
%                     set(gca,'FontSize', FONTSIZE) % Set text size of current axis
%                 end
%                 
%                 % 9-16-22: Plot raw value of wave 1 amplitude change for group
%                 % average
%                 nexttile(t)
%                 n_pts = zeros(n_conditions, 1);
%                 handles = [];
%                 handle_count = 1;
%                 
%                 % Initialize control condition wave 1 amplitude change
%                 % statistics for comparison with noise condition
%                 control_n = zeros(n_timepoints - 1, num_levels);
%                 control_mean = zeros(n_timepoints - 1, num_levels);
%                 control_ste = zeros(n_timepoints - 1, num_levels);
%                 for cc = 1:n_conditions % Iterate over control then noise mice
%                     
%                     hold on
%                     y_pre = [];
%                     for ii = 1:n_timepoints
%                         combined_filter = is_side & is_ingroup & control_noise_selector{cc} & isfreq & metric_selector{mm} & timepoints{ii};
%                         data_selected = finalTable_metadata(combined_filter, :);
%             %             x = data_selected.Properties.VariableNames(COLS_METRICVALS);
%                         x = finalTable_metadata(1, COLS_METRICVALS);
%                         y = data_selected(:, COLS_METRICVALS);
% 
%                         n_pts(fff, ii) = sum(combined_filter);
%                         ymean = nanmean(y{:,:}, 1);
%         %                 yste = nanstd(y{:,:}, 1)/sqrt(n_pts(gg, fff, ii)); % standard error
%                         yste = nanstd(y{:,:}, 1); % Standard deviation
% 
% 
% 
%                         % sort points in order of ascending x value to avoid criss
%                         % crossing of plot lines.
%             %             x_double = cellfun(@str2double, x); % convert cell array of char vectors to double matrix
%                         x_double = x.(1); % convert table 1x1 to double matrix
%                         [sortedX, sortIndex] = sort(x_double);
%                         sortedYmean = ymean(sortIndex);
%                         sortedYstd = yste(sortIndex);
% 
%                         % Remove nan values to plot uninterrupted curve
%                         idx = ~(isnan(sortedYmean));
%                         sortedX_real = sortedX(idx);
%                         sortedYmean_real = sortedYmean(idx);
%                         sortedYstd_real = sortedYstd(idx);
%                         num_levels = numel(sortedX);
%                        
% %                         errorbar(sortedX_real(isdiv10), sortedYmean_real(isdiv10), sortedYstd_real(isdiv10), colors(ii), 'LineWidth', LINEWIDTH)
%         %                 plot(sortedX_real, sortedYmean_real, colors(ii), 'LineWidth', LINEWIDTH)
% 
%                         % Determine if pre and post are significantly different
%                         % at each stimulus level
%                         if ii==1
%                             y_pre = y{:,:};
%                             y_pre = y_pre(:, sortIndex);
%                         else
%                             y_post = y{:,:};
%                             y_post = y_post(:, sortIndex);
%                             d_prime = zeros(num_levels, 1);
%                             delta_wave1amp = zeros(num_levels, 1);
%                             ste_delta_wave1amp = zeros(num_levels, 1);
%                             for zz = 1:num_levels % iterate thru stimulus levels
%                                 this_prey = y_pre(:, zz);
%                                 this_posty = y_post(:, zz);
% 
%                                 if all(isnan(this_prey)) || all(isnan(this_posty))
%                                     continue
%                                 end
% 
%                                 % Remove NaN values
%                                 this_prey = this_prey(~isnan(this_prey));
%                                 this_posty = this_posty(~isnan(this_posty));
%                                 
%                                 % Calculate D' for this stimulus level
%                                 n1 = numel(this_prey);
%                                 mean1 = mean(this_prey);
%                                 std1 = std(this_prey);
%                                 n2 = numel(this_posty);
%                                 mean2 = mean(this_posty);
%                                 std2 = std(this_posty);
%                                 [npool,meanpool,stdpool] = pooledmeanstd(n1,mean1,std1,n2,mean2,std2);
%                                 d_prime(zz) = (mean2 - mean1)/stdpool;
%                                 
%                                 % 9-16-22 Also calculate change in wave 1
%                                 % amplitude compared to baseline
%                                 delta_wave1amp(zz) = mean2 - mean1;
%                                 ste1 = std1/sqrt(n1);
%                                 ste2 = std2/sqrt(n2);
%                                 ste_delta_wave1amp(zz) = sqrt(ste1^2 + ste2^2);
% 
%                                 % Compare wave 1 amplitude change between control and noise conditions using 2-sample unpaired t-test
%                                 if strcmpi(CONTROL_NOISE_LABEL{cc}, 'control') % control condition - store wave 1 amplitude changes
%                                     control_n(ii - 1, zz) = n2;
%                                     control_mean(ii - 1, zz) = delta_wave1amp(zz);
%                                     control_ste(ii - 1, zz) = ste_delta_wave1amp(zz);
%                                 else  % noise condition - perform unpaired 2 sample t-test to determine if different from control condition
%                                     is_different = ttest2_mean_ste(control_n(ii - 1, zz), control_mean(ii - 1, zz),...
%                                         control_ste(ii - 1, zz), n2, delta_wave1amp(zz), ste_delta_wave1amp(zz));
%                                     if is_different < P_CRIT && mod(sortedX(zz), 10)==0
%                                         % asterisk on noise condition's
%                                         % delta wave 1 amp curve if
%                                         % significantly different from
%                                         % control change
%                                         text(sortedX(zz), delta_wave1amp(zz), '*', ...
%                                             'Color', colors(ii), 'FontSize', FONTSIZE*2, 'FontWeight', 'bold') 
%                                     end
%                                 end
%                             end
%                             
%                             % Remove stimulus levels that are not divisible by 10
%                             isdiv10 = mod(sortedX, 10)==0;
%                             
%                             % Plot raw wave 1 amplitude change for this post-exposure date
%                             if strcmpi(CONTROL_NOISE_LABEL{cc}, 'control') % control lines are dashed
%                                 handles(handle_count) = errorbar(sortedX(isdiv10), delta_wave1amp(isdiv10), ...
%                                     ste_delta_wave1amp(isdiv10), colors(ii), 'LineWidth', LINEWIDTH, 'LineStyle', "--");
%                             else % noise exposed mouse lines are solid
%                                 handles(handle_count) = errorbar(sortedX(isdiv10), delta_wave1amp(isdiv10), ...
%                                     ste_delta_wave1amp(isdiv10), colors(ii), 'LineWidth', LINEWIDTH);
%                             end
%                             
%                             handle_count = handle_count + 1;
%                         end
% 
%                     end
% 
%                     %TODO: add statistical significance test between pre and post
%                     %noise curves
% 
%                     hold off
%                     legend(handles, [strcat(NAME_TIMEPOINTS(2:end), '-control'), strcat(NAME_TIMEPOINTS(2:end), '-noise')], 'Location','northwest', 'AutoUpdate','off')
%                     yline(0, ':')
%                     yline(-1, ':')
%                     yline(1, ':')
% 
%                     xlabel('Amplitude (dB)')
%                     ylabel(YLABEL_METRIC_RAWCHANGE{mm})
%                     ylim(YLIM_METRIC_RAW{mm})
%         %             ylim([-0.1, 1.1]) % For normalized wave1amp
%                     xlim([-30, 100])
%             %         title([METRIC_NAMES{mm}, ' ', ', ', num2str(FREQ(fff)), ' Hz. N=', num2str(n_pts(fff, 1))]); % TODO: fix n_pts, which is currently 33 times (n_timepts*num_levels) the actual number of mice
% %                     title([CONTROL_NOISE_LABEL{cc}, ', ', num2str(FREQ(fff)), ' Hz. N=', num2str(n_pts(fff, ii)), ' ears'], 'Interpreter','none');
%                     set(gca,'FontSize', FONTSIZE) % Set text size of current axis
%                 end
% 
%             % Label the tiled layout figure
%     %         t.TileSpacing = 'none';
%             t.TileSpacing = 'tight';
%             t.Padding = 'tight';
%             plot_title = [GROUP_NAMES{gg}, ', ', SIDES{nn}, ' side, ', METRIC_NAMES{mm}];
%             title(t, plot_title, 'FontSize', 24)
%     %         ylabel(t, 'Amplitude (nV)', 'FontSize', 24)
%     %         xlabel(t, 'Time (ms)', 'FontSize', 24)
% 
%             % Maximize figure window size
%             fig.WindowState = 'maximized';
% 
%             % Save figure to results folder
%             save_filename = ['levelfunction_group', GROUP_NAMES{gg}, '_', SIDES{nn}, '_', METRIC_NAMES{mm}];
%             savefig_multiformat(gcf, SAVE_PATH, save_filename)
%         end
%     end
% end


disp('Done.')