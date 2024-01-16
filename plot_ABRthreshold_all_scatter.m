% plot_ABRthreshold_all_scatter.m
%
% Create scatter plot of calculated ABR thresholds that were output in
% 2021MMDDDD_abr_output.xlsx files from wrapper_analyzeABR_all.m.
%
% Adapted from plot_ABRthreshold_levelfunctions.m
%
% 11/9/21 George Liu
% Dependencies: load_all_ABR_DPOAE_analysis.m, get_mousefile_metadata.m, savefig_multiformat.m

%% Constants
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures
LOAD_PATH_ABR = 'd:\users\admin\Documents\George\Results_12-7-21_abr';
LOAD_PATH_DPOAE = 'd:\users\admin\Documents\George\Results_11-9-21_dpoae';
% BATCHGROUPS = {1:4, 5:6, 7, 8, 9, 1:9}; % group sequentially by date of experimetal batches
% GROUP_NAMES = {'1', '2', '3', '4', '5', '6'};
BATCHGROUPS = {5:6, [2:4,7,9], 8, 2:9}; % group by noise level: 96 dB, 97 dB, and 99 dB respectively.
% BATCHGROUPS = {6, [2:4,7,9], 8, [2:4,6:9]}; % group by noise level: 96 dB, 97 dB, and 99 dB respectively. Remove group 5 (L and R ears asymmetric).
GROUP_NAMES = {'96 dB', '97 dB', '99 dB', 'all'};
% BATCHGROUPS = {2:9}; % group by noise level: 96 dB, 97 dB, and 99 dB respectively.
% GROUP_NAMES = {'all'};
N_GROUPS = numel(BATCHGROUPS);
% METRIC_NAMES = {'Liberman', 'Oghalai', 'Innerprod', 'Innerprod_AUC', 'D_rms', 'D_z'};
METRIC_NAMES = {'Liberman', 'Oghalai', 'Innerprod', 'Innerprod_AUC'};
n_metrics = numel(METRIC_NAMES);
FREQ = [8000, 16000, 32000];
n_freq = numel(FREQ);
COL_METRICVALS_FIRSTCOL = 24; % first column index of metric values columns in finalTable_metadata table
CONTROL_NOISE_LABEL = {'control', 'noise'};
YLIM_METRIC = {[-0.4, 1], [0, 6250], [0, 2.5*10^-9], [0.5, 1]};
YLABEL_METRIC = {'Correlation coefficient (a.u.)', 'Max p-p amplitude (nV)', 'Mean inner product (nv^2)', 'AUC (a.u.)'}; 
SIDES = ["left", "right"];
n_sides = numel(SIDES);
SIDES_SYMBOLS = ['x', 'o'];
NAME_TIMEPOINTS = ["Pre", "Post 24h", "Post 1w"];
P_CRIT = 0.05;
LINEWIDTH = 2;
FONTSIZE = 14;

% METRIC_NAMES = {'Wave1amp'};
% n_metrics = numel(METRIC_NAMES);
% YLIM_METRIC = {[0, 6250]};
% YLABEL_METRIC = {'Wave 1 amplitude (nV)'}; 
% 
% METRIC_NAMES = {'Innerprod_window'};
% n_metrics = numel(METRIC_NAMES);
% YLIM_METRIC = {[-1*10^-12, 1*10^-11]};
% YLABEL_METRIC = {'Inner prod window (nV^2)'}; 


%% Load data from all ABR analysis output table files into one aggregate table
is_abr = true;
finalTable = load_all_ABR_DPOAE_analysis(LOAD_PATH_ABR, is_abr);


%% Add metadata to table
n_rows = size(finalTable, 1);
for k=1:n_rows
    disp(['Working on metadata row ', num2str(k), ' out of ', num2str(n_rows)])
    date_name = finalTable.Filenames(k);
    [date, name, studytype, side, metadata] = get_mousefile_metadata(date_name);
    metadata.('Date') = date;
    metadata.('Side') = side;
    metadata.('Studytype') = studytype;
    if k==1
        metadata_aggregate = metadata;
    else
        metadata_aggregate = vertcat(metadata_aggregate,metadata);
    end
end

finalTable_metadata = [metadata_aggregate, finalTable];

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

%% Plot ABR threshold shift per frequency for noise exposed mice. 
% Plot all threshold methods for all groups of mice.
colors = ['b', 'g', 'r'];

for nn = 1:numel(CONTROL_NOISE_LABEL) 
    % Iterate over control then noise group

    for mm = 1:n_metrics
        % get metric
        is_metric = metric_selector{mm};
        
        % New tiled plot for each metric.
        fig = figure;
        ht = tiledlayout(3, N_GROUPS, 'TileIndexing', 'columnmajor');
%             ht = tiledlayout(1, 2, 'TileIndexing', 'columnmajor');
        
        for gg = 1:N_GROUPS
            %get group
            is_group = any(finalTable_metadata.('Batch') == BATCHGROUPS{gg}, 2);

            % Plot thresholds vs freq at each time point
            nexttile(ht)
            hold on
            
            x_all = cell(1, n_timepoints);
            y_all = cell(1, n_timepoints);
            ymean_all = zeros(n_freq, n_timepoints);
            ystd = zeros(n_freq, n_timepoints);
            y_rel_std = zeros(n_freq, n_timepoints);
            p = zeros(1, n_timepoints);
            stats_freq = cell(1, n_timepoints);
            h_multcomp = cell(1, n_timepoints);
            n_pts = zeros(n_freq, n_timepoints);
            for ii = 1:n_timepoints
                data_selected = finalTable_metadata(is_group & is_metric & control_noise_selector{nn} & timepoints{ii}, :);
                x = data_selected.('Frequency');
                y = data_selected.('Threshold');

                ymean = zeros(n_freq, 1);
                yste = zeros(n_freq, 1);
                for ff = 1:n_freq
                    n_pts(ff, ii) = sum(x==FREQ(ff));
                    ymean(ff) = nanmean(y(x==FREQ(ff)));
                    ystd(ff, ii) = nanstd(y(x==FREQ(ff)));
                    yste(ff) = ystd(ff, ii)/sqrt(n_pts(ff, ii));
                    y_rel_std(ff, ii) = ystd(ff)/abs(ymean(ff));
                end
                
                % Plot
                errorbar(FREQ, ymean, yste, colors(ii), 'LineWidth', LINEWIDTH)
                
%                 % Plot asterisk next to plot if 1 way anova shows that
%                 % threshold is different across frequencies
%                 [p(ii), tbl, stats_freq{ii}] = anova1(y,x, 'off');
%                 if p(ii) < P_CRIT
%                     text(34000, ymean(ff), '#', 'Color', colors(ii), 'FontSize', FONTSIZE, 'FontWeight', 'bold')
%                 end
                
%                 [c, m, h_multcomp{ii}] = multcompare(stats_freq{ii}); % determine which frequency is outlier

                % Cache mean threshold values to later analyze thresholds
                % vs time
                x_all{ii} = x;
                y_all{ii} = y;
                ymean_all(:, ii) = ymean;
            end
            
            % Check threshold stability vs time
            for ff = 1:n_freq
                thresh_ntimepoints = cell(1, n_timepoints);
                for ii = 1:n_timepoints
                    this_x = x_all{ii};
                    this_y = y_all{ii};
                    thresh_ntimepoints{ii} = this_y(this_x==FREQ(ff));
                end
                thresh_ntimepoints_vec = [thresh_ntimepoints{1}; thresh_ntimepoints{2}; thresh_ntimepoints{3}];
                thresh_ntimepoints_group = [repmat(NAME_TIMEPOINTS(1), numel(thresh_ntimepoints{1}), 1); ...
                    repmat(NAME_TIMEPOINTS(2), numel(thresh_ntimepoints{2}), 1); ...
                    repmat(NAME_TIMEPOINTS(3), numel(thresh_ntimepoints{3}), 1)];
                p2 = anova1(thresh_ntimepoints_vec, thresh_ntimepoints_group, 'off');
                
                if p2 < P_CRIT
                    [max_ymean, ind] = max(ymean_all(ff, :));
                    text(FREQ(ff), max_ymean + 5, '*', 'Color', colors(ind), 'FontSize', FONTSIZE, 'FontWeight', 'bold') % asterisk above highest point at frequency if ANOVA positive
                end
                
                [~, p3] = ttest2(thresh_ntimepoints{1}, thresh_ntimepoints{3}); % compare pre and post 1 w threshold group means
                if p3 < P_CRIT && ff==n_freq
                    text(33000, ymean(ff), '*', 'Color', colors(ii), 'FontSize', FONTSIZE, 'FontWeight', 'bold')
                end
                    
            end
            
            hold off
            if gg==1
                legend(NAME_TIMEPOINTS, 'Location','northwest')
            else
                set(gca,'yticklabel',[]) % Remove y-axis labels and legend from all but first column.
            end
            xlabel('Frequency (Hz)')
%             ylabel('Threshold (dB)')
            ylim([-10, 90])
            xlim([5000, 35000])
            set(gca,'FontSize', FONTSIZE) % Set text size of current axis
            title([METRIC_NAMES{mm}, ' ', CONTROL_NOISE_LABEL{nn}, ', ',  GROUP_NAMES{gg}, '. N=', num2str(n_pts(1)), ' ears']);
    %         disp(n_pts)
            
            % Plot standard deviation
            nexttile
            hold on
            for ii = 1:n_timepoints
                plot(FREQ, ystd(:, ii), colors(ii), 'LineWidth', LINEWIDTH)
            end
            hold off
            legend(NAME_TIMEPOINTS, 'Location','northwest')
            xlabel('Frequency (Hz)')
            ylabel('Standard deviation (dB)')
            ylim([0, 90])
            title([METRIC_NAMES{mm}, ' ', CONTROL_NOISE_LABEL{nn}, ', ', GROUP_NAMES{gg}, '. N=', num2str(n_pts(1)), ' ears']);
            set(gca,'FontSize', FONTSIZE) % Set text size of current axis
            
            %Plot relative standard deviation 
            nexttile
            hold on
            for ii = 1:n_timepoints
                plot(FREQ, y_rel_std(:, ii), colors(ii), 'LineWidth', LINEWIDTH)
            end
            hold off
            legend(NAME_TIMEPOINTS, 'Location','northwest')
            xlabel('Frequency (Hz)')
            ylabel('Coefficient of variation (a.u.)')
            ylim([0, 5])
            title([METRIC_NAMES{mm}, ' ', CONTROL_NOISE_LABEL{nn}, ', ', GROUP_NAMES{gg}, '. N=', num2str(n_pts(1)), ' ears']);
            set(gca,'FontSize', FONTSIZE) % Set text size of current axis
            
%             % Plot ANOVA 1 multiple comparison results for frequency
%             nexttile
            
            % Final tiled chart formatting
            ht.TileSpacing = 'tight';
            ht.Padding = 'tight';
%             ht.TileSpacing = 'tight';
%             ht.Padding = 'none';
            % Maximize figure window size
            fig.WindowState = 'maximized';
%             set(gcf, 'Position', [434, 385, 1191, 420])

        end
        
        % Save figure to results folder
        save_filename = ['threshold_shift_', METRIC_NAMES{mm}, '_', CONTROL_NOISE_LABEL{nn}];
        savefig_multiformat(gcf, SAVE_PATH, save_filename)
    end
end

%% 12-7-21: Plot ABR threshold shift per individual mice (per frequency). 
% Plot change in threshold from pre to post-24h and post-1w for all groups of mice.
colors = ['b', 'g', 'r'];

for nn = 1:numel(CONTROL_NOISE_LABEL) 
    % Iterate over control then noise group

    for mm = 1:n_metrics
        % get metric
        is_metric = metric_selector{mm};
        
        % New tiled plot for each metric.
        fig = figure;
        ht = tiledlayout(3, N_GROUPS, 'TileIndexing', 'columnmajor');
%             ht = tiledlayout(1, 2, 'TileIndexing', 'columnmajor');
        
        for gg = 1:N_GROUPS
            %get group
            is_group = any(finalTable_metadata.('Batch') == BATCHGROUPS{gg}, 2);

            data_selected = finalTable_metadata(is_group & is_metric & control_noise_selector{nn}, :);
            number_mice = finalTable_metadata.('Number');
            unique_numbers_mice = unique(number_mice);
            n_numbers = length(unique_numbers_mice);
            
            thresh = NaN(n_freq, n_numbers, n_sides, n_timepoints);
            thresh_shift = NaN(n_freq, n_numbers, n_sides, n_timepoints);
            

            % Organize threshold data by mice
            for ff = 1:n_freq
                is_freq = finalTable_metadata.('Frequency') == FREQ(ff);
                for zz = 1:n_numbers
                    is_number = number_mice == unique_numbers_mice(zz);
                    for ss = 1:n_sides
                        is_side = strcmpi(finalTable_metadata.('Side'), SIDES(ss));
                        for ii = 1:n_timepoints
                            data_selected_freq_number = finalTable_metadata(is_group & is_metric & control_noise_selector{nn} & is_freq & is_number & is_side & timepoints{ii}, :);
                            if ~isempty(data_selected_freq_number)
                                thresh(ff, zz, ss, ii) = data_selected_freq_number.('Threshold'); 
                            end

                            % Calculate difference of post thresholds from baseline
                            % threshold
                            thresh_shift(ff, zz, ss, ii) = thresh(ff, zz, ss, ii) - thresh(ff, zz, ss, 1);
                        end
                    end
                end
            end
            
            % Plot thresholds shifts vs freq for post-exposure times
            nexttile(ht)
            hold on
            
            x_all = cell(1, n_timepoints);
            y_all = cell(1, n_timepoints);
            ymean_all = zeros(n_freq, n_timepoints);
            ystd = zeros(n_freq, n_timepoints);
            y_rel_std = zeros(n_freq, n_timepoints);
            p = zeros(1, n_timepoints);
            stats_freq = cell(1, n_timepoints);
            h_multcomp = cell(1, n_timepoints);
            n_pts = zeros(n_freq, n_timepoints);
            
%             ymean = zeros(n_freq, 1);
%             ystd = zeros(n_freq, 1);
            yste = zeros(n_freq, 1);
            for ii = 1:n_timepoints
                x = FREQ;
                ymean = nanmean(thresh_shift(:, :, :, ii), [2, 3]);
                
                for ff = 1:n_freq
                    this_threshshift = thresh_shift(ff, :, :, ii);
                    ystd(ff, ii) = nanstd(this_threshshift(:));
                    n_pts(ff, ii) = sum(~isnan(this_threshshift(:)));
                    yste(ff) = ystd(ff, ii)./sqrt(n_pts(ff, ii));
                    y_rel_std(ff, ii) = ystd(ff, ii) ./ ymean(ff);
                end

                % Plot
                errorbar(FREQ, ymean, yste, colors(ii), 'LineWidth', LINEWIDTH)
                
%                 % Plot asterisk next to plot if 1 way anova shows that
%                 % threshold is different across frequencies
%                 [p(ii), tbl, stats_freq{ii}] = anova1(y,x, 'off');
%                 if p(ii) < P_CRIT
%                     text(34000, ymean(ff), '#', 'Color', colors(ii), 'FontSize', FONTSIZE, 'FontWeight', 'bold')
%                 end
                
%                 [c, m, h_multcomp{ii}] = multcompare(stats_freq{ii}); % determine which frequency is outlier

                % Cache mean threshold values to later analyze thresholds
                % vs time
                x_all{ii} = x;
                y_all{ii} = thresh_shift(:, :, :, ii);
                ymean_all(:, ii) = ymean;
            end
            
%             % Check threshold stability vs time
%             for ff = 1:n_freq
%                 thresh_ntimepoints = cell(1, n_timepoints);
%                 for ii = 1:n_timepoints
%                     this_x = x_all{ii};
%                     this_y = y_all{ii};
%                     thresh_ntimepoints{ii} = this_y(this_x==FREQ(ff));
%                 end
%                 thresh_ntimepoints_vec = [thresh_ntimepoints{1}; thresh_ntimepoints{2}; thresh_ntimepoints{3}];
%                 thresh_ntimepoints_group = [repmat(NAME_TIMEPOINTS(1), numel(thresh_ntimepoints{1}), 1); ...
%                     repmat(NAME_TIMEPOINTS(2), numel(thresh_ntimepoints{2}), 1); ...
%                     repmat(NAME_TIMEPOINTS(3), numel(thresh_ntimepoints{3}), 1)];
%                 p2 = anova1(thresh_ntimepoints_vec, thresh_ntimepoints_group, 'off');
%                 
%                 if p2 < P_CRIT
%                     [max_ymean, ind] = max(ymean_all(ff, :));
%                     text(FREQ(ff), max_ymean + 5, '*', 'Color', colors(ind), 'FontSize', FONTSIZE, 'FontWeight', 'bold') % asterisk above highest point at frequency if ANOVA positive
%                 end
%                 
%                 [~, p3] = ttest2(thresh_ntimepoints{1}, thresh_ntimepoints{3}); % compare pre and post 1 w threshold group means
%                 if p3 < P_CRIT && ff==n_freq
%                     text(33000, ymean(ff), '*', 'Color', colors(ii), 'FontSize', FONTSIZE, 'FontWeight', 'bold')
%                 end
%                     
%             end
            
            hold off
            if gg==1
                legend(NAME_TIMEPOINTS, 'Location','northwest')
            else
                set(gca,'yticklabel',[]) % Remove y-axis labels and legend from all but first column.
            end
            xlabel('Frequency (Hz)')
%             ylabel('Threshold (dB)')
            ylim([-10, 90])
            xlim([5000, 35000])
            set(gca,'FontSize', FONTSIZE) % Set text size of current axis
            title([METRIC_NAMES{mm}, ' ', CONTROL_NOISE_LABEL{nn}, ', ',  GROUP_NAMES{gg}, '. N=', num2str(n_pts(1)), ' ears']);
    %         disp(n_pts)
            
            % Plot standard deviation
            nexttile
            hold on
            for ii = 1:n_timepoints
                plot(FREQ, ystd(:, ii), colors(ii), 'LineWidth', LINEWIDTH)
            end
            hold off
            legend(NAME_TIMEPOINTS, 'Location','northwest')
            xlabel('Frequency (Hz)')
            ylabel('Standard deviation (dB)')
            ylim([0, 90])
            title([METRIC_NAMES{mm}, ' ', CONTROL_NOISE_LABEL{nn}, ', ', GROUP_NAMES{gg}, '. N=', num2str(n_pts(1)), ' ears']);
            set(gca,'FontSize', FONTSIZE) % Set text size of current axis
            
            %Plot relative standard deviation 
            nexttile
            hold on
            for ii = 1:n_timepoints
                plot(FREQ, y_rel_std(:, ii), colors(ii), 'LineWidth', LINEWIDTH)
            end
            hold off
            legend(NAME_TIMEPOINTS, 'Location','northwest')
            xlabel('Frequency (Hz)')
            ylabel('Coefficient of variation (a.u.)')
            ylim([0, 10])
            title([METRIC_NAMES{mm}, ' ', CONTROL_NOISE_LABEL{nn}, ', ', GROUP_NAMES{gg}, '. N=', num2str(n_pts(1)), ' ears']);
            set(gca,'FontSize', FONTSIZE) % Set text size of current axis
            
%             % Plot ANOVA 1 multiple comparison results for frequency
%             nexttile
            
            % Final tiled chart formatting
            ht.TileSpacing = 'tight';
            ht.Padding = 'tight';
%             ht.TileSpacing = 'tight';
%             ht.Padding = 'none';
            % Maximize figure window size
            fig.WindowState = 'maximized';
%             set(gcf, 'Position', [434, 385, 1191, 420])

        end
        
        % Save figure to results folder
        save_filename = ['threshold_shift_', METRIC_NAMES{mm}, '_', CONTROL_NOISE_LABEL{nn}];
        savefig_multiformat(gcf, SAVE_PATH, save_filename)
    end
end



%% Plot left ear ABR threshold vs right ear ABR threshold
is_left = strcmpi(finalTable_metadata.('Side'), SIDES(1));

% Iterate over control then noise group
for nn = 1:numel(CONTROL_NOISE_LABEL) 
    
    for gg = 1:N_GROUPS %get group   
        is_group = any(finalTable_metadata.('Batch') == BATCHGROUPS{gg}, 2);
        
        fig = figure;
        ht = tiledlayout(2, n_metrics, 'TileIndexing', 'columnmajor');

        for mm = 1:n_metrics
            % get metric
            is_metric = metric_selector{mm};
            data_selected_left = finalTable_metadata(is_metric & control_noise_selector{nn} & is_group & is_left, :);
            data_selected_right = finalTable_metadata(is_metric & control_noise_selector{nn} & is_group & ~is_left, :);

            % create list of ABR thresholds for which left and right ear
            % measurements exist
            thresh_left = [];
            thresh_right = [];
            for z = 1:size(data_selected_left, 1)
                % Check that ABR threshold for same mouse, frequency, and date also exists for right ear
                this_mouse_name = data_selected_left.('Mouse_name'){z};
                this_mouse_date = data_selected_left.('Date'){z};
                this_freq = data_selected_left.('Frequency')(z);
                is_right_mouse = strcmpi(data_selected_right.('Mouse_name'), this_mouse_name);
                is_right_date = strcmpi(data_selected_right.('Date'), this_mouse_date);
                is_right_freq = data_selected_right.('Frequency')==this_freq;
                is_matching_row_right = data_selected_right(is_right_mouse & is_right_date & is_right_freq, :);
                if size(is_matching_row_right, 1)==1
                    thresh_left = [thresh_left, data_selected_left.('Threshold')(z)];
                    thresh_right = [thresh_right, is_matching_row_right.('Threshold')(1)];
                end
            end
            nexttile
            if mm>2 % add random jitter to innerproduct thresholds so points don't all overlap
                rand_jitter = 2*(rand(size(thresh_left))-0.5);
                s = scatter(thresh_right + rand_jitter, thresh_left + rand_jitter, 'filled');
            else
                s = scatter(thresh_right, thresh_left, 'filled');
            end
            s.AlphaData = 0.1*ones(size(thresh_right)); % set transparency of filled circle marks
            s.MarkerFaceAlpha = 'flat';
            % plot line y=x
            x=-90:90;
            hold on, plot(x, x, '--', 'LineWidth', LINEWIDTH), hold off

            axis equal
            ylim([-40, 90])
            xlim([-40, 90])
            xlabel('Right ear ABR threshold (dB)')
            title(METRIC_NAMES{mm})
            set(gca,'FontSize', FONTSIZE) % Set text size of current axis
            if mm==1
                ylabel('Left ear ABR threshold (dB)')
%             else
%                 set(gca,'yticklabel',[]) % Remove y-axis labels and legend from all but first column.
            end

            % Plot Bland-Altman plot
            left_minus_right_thresh = thresh_left - thresh_right;
            ave_thresh = (thresh_left + thresh_right)/2;
            nexttile
            if mm>2
                s=scatter(ave_thresh + rand_jitter, left_minus_right_thresh, 'filled');
            else
                s=scatter(ave_thresh, left_minus_right_thresh, 'filled');
            end
            s.AlphaData = 0.1*ones(size(ave_thresh)); % set transparency of filled circle marks
            s.MarkerFaceAlpha = 'flat';
            xlim([-40, 90])
            ylim([-100, 100])
            limits_of_agreement = [-1.96*nanstd(left_minus_right_thresh), 1.96*nanstd(left_minus_right_thresh)];
            hold on
            yline(limits_of_agreement)
            hold off
            xlabel('Average of left and right thresholds')
            if mm==1
                ylabel('Left minus right threshold')
%             else
%                 set(gca,'yticklabel',[]) % Remove y-axis labels and legend from all but first column.
            end

            title(ht, [CONTROL_NOISE_LABEL{nn}, ', group ', GROUP_NAMES{gg}])
            set(gca,'FontSize', FONTSIZE) % Set text size of current axis
            
            % Final tiled chart formatting
            ht.TileSpacing = 'tight';
            ht.Padding = 'tight';
            % Maximize figure window size
            fig.WindowState = 'maximized';
            
            % Save figure to results folder
            save_filename = ['left_right_', CONTROL_NOISE_LABEL{nn}, '_group', GROUP_NAMES{gg}];
            savefig_multiformat(gcf, SAVE_PATH, save_filename)
        end
    end
end

%% Plot scatter of threshold vs date of ABR experiment
NAME_TIMEPOINTS_SIDES = ["L Pre", "L 24h", "L 1w", "R Pre", "R 24h", "R 1w"];
                    
for nn = 1:numel(CONTROL_NOISE_LABEL)
    fig = figure;
    ht = tiledlayout(n_freq, n_metrics-1, 'TileIndexing', 'columnmajor');
    for mm = 1:n_metrics-1
        is_metric = metric_selector{mm};
        this_filter = is_metric & control_noise_selector{nn};
        data_selected = finalTable_metadata(this_filter, :);
        z = data_selected.('Frequency');
        y = data_selected.('Threshold');
        x = data_selected.('Date');
        t = datetime(x, 'InputFormat', 'yyyyMMdd');

        for ff = 1:n_freq
            nexttile(ht)
            hold on
            for ss = 1:n_sides
                for ii = 1:n_timepoints
                    is_timepoint = timepoints{ii};
                    is_freq = z==FREQ(ff);
                    is_side = strcmp(data_selected.('Side'), SIDES(ss));
                    second_filter = is_freq & is_timepoint(this_filter) & is_side;

                    t_subset = t(second_filter);
                    y_subset = y(second_filter);
                    scatter(t_subset, y_subset, [], colors(ii), SIDES_SYMBOLS(ss))
                end
            end
            hold off
            if ff==1 && mm==1
                legend(NAME_TIMEPOINTS_SIDES, 'Location', 'northeast', 'NumColumns', 2)
            end
            xlabel('Date')
            ylabel('Threshold (dB)')
            ylim([-30, 90])
            title([METRIC_NAMES{mm}, ' ', CONTROL_NOISE_LABEL{nn}, ' ', num2str(FREQ(ff)), ' Hz'])
            set(gca,'FontSize', FONTSIZE) % Set text size of current axis
        end
    end

    ht.TileSpacing = 'tight';
    ht.Padding = 'tight';
    %     plot_title = ['ABR @ ', num2str(this_freq), ' Hz'];
    %     title(t, plot_title)
    %     ylabel(t, 'Amplitude (nV)')

    xlabel(ht, 'Date', 'FontSize', FONTSIZE)
    ylabel(ht, 'Threshold (dB)', 'FontSize', FONTSIZE)
    % Maximize figure window size
    fig.WindowState = 'maximized';
    
    % Save figure to results folder
    save_filename = ['date_vs_threshold_', METRIC_NAMES{mm}, '_', CONTROL_NOISE_LABEL{nn}, '_', GROUP_NAMES{gg}];
    savefig_multiformat(gcf, SAVE_PATH, save_filename)
end

%% Load analyzed DPOAE thresholds
is_abr = false;
finalTable_dpoae = load_all_ABR_DPOAE_analysis(LOAD_PATH_DPOAE, is_abr);

% Add metadata to DPOAE output table
n_rows = size(finalTable_dpoae, 1);
for k=1:n_rows
    disp(['Working on metadata row ', num2str(k), ' out of ', num2str(n_rows)])
    date_name = finalTable_dpoae.Filenames(k);
    [date, name, studytype, side, metadata] = get_mousefile_metadata(date_name);
    metadata.('Date') = date;
    metadata.('Side') = side;
    metadata.('Studytype') = studytype;
    if k==1
        metadata_aggregate_dpoae = metadata;
    else
        metadata_aggregate_dpoae = vertcat(metadata_aggregate_dpoae, metadata);
    end
end

finalTable_dpoae_metadata = [metadata_aggregate_dpoae, finalTable_dpoae];

% Create filters
% Noise exposure
is_noise_dpoae = logical(finalTable_dpoae_metadata.('Is_noise_exposed'));
is_control_dpoae = ~is_noise_dpoae;
control_noise_selector_dpoae = {is_control_dpoae, is_noise_dpoae};

% Timepoint
is_pre_dpoae = strcmp(string(finalTable_dpoae_metadata.('Date_1')), string(finalTable_dpoae_metadata.('Date')));
is_post24_dpoae = strcmp(string(finalTable_dpoae_metadata.('Date_2')), string(finalTable_dpoae_metadata.('Date')));
is_post1w_dpoae = strcmp(string(finalTable_dpoae_metadata.('Date_3')), string(finalTable_dpoae_metadata.('Date')));
timepoints_dpoae = {is_pre_dpoae, is_post24_dpoae, is_post1w_dpoae};
n_timepoints_dpoae = numel(timepoints_dpoae);

%% Plot DPOAE threshold shift per frequency for noise exposed mice. 
% Plot all threshold methods for all groups of mice.
colors = ['b', 'g', 'r'];

for nn = 1:numel(CONTROL_NOISE_LABEL) 
    % Iterate over control then noise group

    % New tiled plot for each metric.
    fig = figure;
    ht = tiledlayout(1, N_GROUPS, 'TileIndexing', 'columnmajor');

    for gg = 1:N_GROUPS
        %get group
        is_group = any(finalTable_dpoae_metadata.('Batch') == BATCHGROUPS{gg}, 2);

        % Plot thresholds vs freq at each time point
        nexttile(ht)
        hold on

        x_all = cell(1, n_timepoints);
        y_all = cell(1, n_timepoints);
        ymean_all = zeros(n_freq, n_timepoints);
        ystd = zeros(n_freq, n_timepoints);
        y_rel_std = zeros(n_freq, n_timepoints);
        p = zeros(1, n_timepoints);
        stats_freq = cell(1, n_timepoints);
        h_multcomp = cell(1, n_timepoints);
        n_pts = zeros(n_freq, n_timepoints);
        for ii = 1:n_timepoints_dpoae
            data_selected = finalTable_dpoae_metadata(is_group & control_noise_selector_dpoae{nn} & timepoints_dpoae{ii}, :);
            x = data_selected.('Frequency');
            y = data_selected.('Threshold');

            ymean = zeros(n_freq, 1);
            yste = zeros(n_freq, 1);
            for ff = 1:n_freq
                n_pts(ff, ii) = sum(x==FREQ(ff));
                ymean(ff) = nanmean(y(x==FREQ(ff)));
                ystd(ff, ii) = nanstd(y(x==FREQ(ff)));
                yste(ff) = ystd(ff, ii)/sqrt(n_pts(ff, ii));
                y_rel_std(ff, ii) = ystd(ff)/abs(ymean(ff));
            end

            % Plot
            errorbar(FREQ, ymean, yste, colors(ii), 'LineWidth', LINEWIDTH)

            % Cache mean threshold values to later analyze thresholds
            % vs time
            x_all{ii} = x;
            y_all{ii} = y;
            ymean_all(:, ii) = ymean;
        end

        % Check threshold stability vs time
        for ff = 1:n_freq
            thresh_ntimepoints = cell(1, n_timepoints);
            for ii = 1:n_timepoints
                this_x = x_all{ii};
                this_y = y_all{ii};
                thresh_ntimepoints{ii} = this_y(this_x==FREQ(ff));
            end
            thresh_ntimepoints_vec = [thresh_ntimepoints{1}; thresh_ntimepoints{2}; thresh_ntimepoints{3}];
            thresh_ntimepoints_group = [repmat(NAME_TIMEPOINTS(1), numel(thresh_ntimepoints{1}), 1); ...
                repmat(NAME_TIMEPOINTS(2), numel(thresh_ntimepoints{2}), 1); ...
                repmat(NAME_TIMEPOINTS(3), numel(thresh_ntimepoints{3}), 1)];
            p2 = anova1(thresh_ntimepoints_vec, thresh_ntimepoints_group, 'off');

            if p2 < P_CRIT
                [max_ymean, ind] = max(ymean_all(ff, :));
                text(FREQ(ff), max_ymean + 5, '*', 'Color', colors(ind), 'FontSize', FONTSIZE, 'FontWeight', 'bold') % asterisk above highest point at frequency if ANOVA positive
            end

            [~, p3] = ttest2(thresh_ntimepoints{1}, thresh_ntimepoints{3}); % compare pre and post 1 w threshold group means
            if p3 < P_CRIT && ff==n_freq
                text(33000, ymean(ff), '*', 'Color', colors(ii), 'FontSize', FONTSIZE, 'FontWeight', 'bold')
            end

        end

        hold off
        if gg==1 
            if nn==1
                legend(NAME_TIMEPOINTS, 'Location','northeast')
            end
            ylabel('Threshold (dB)')
        else
            set(gca,'yticklabel',[]) % Remove y-axis labels and legend from all but first column.
        end
        xlabel('Frequency (Hz)')
        ylim([30, 90])
        xlim([5000, 35000])
        set(gca,'FontSize', FONTSIZE) % Set text size of current axis
        title([CONTROL_NOISE_LABEL{nn}, ', ',  GROUP_NAMES{gg}, '. N=', num2str(n_pts(1)), ' ears']);
    end
    
    % Final tiled chart formatting
    ht.TileSpacing = 'tight';
    ht.Padding = 'tight';
    set(gcf, 'Position', [434, 385, 1191, 420])

    % Save figure to results folder
    save_filename = ['DPOAEthreshold_shift_', CONTROL_NOISE_LABEL{nn}];
    savefig_multiformat(gcf, SAVE_PATH, save_filename)
end

%% Plot scatter of DPOAE threshold vs date 
NAME_TIMEPOINTS_SIDES = ["L Pre", "L 24h", "L 1w", "R Pre", "R 24h", "R 1w"];
                    
for nn = 1:numel(CONTROL_NOISE_LABEL)
    fig = figure;
    ht = tiledlayout(1, n_freq, 'TileIndexing', 'columnmajor');
    noise_filter = control_noise_selector_dpoae{nn};
    data_selected = finalTable_dpoae_metadata(noise_filter, :);
    z = data_selected.('Frequency');
    y = data_selected.('Threshold');
    x = data_selected.('Date');
    t = datetime(x, 'InputFormat', 'yyyyMMdd');

    for ff = 1:n_freq
        nexttile(ht)
        hold on
        for ss = 1:n_sides
            for ii = 1:n_timepoints
                is_timepoint = timepoints_dpoae{ii};
                is_freq = z==FREQ(ff);
                is_side = strcmp(data_selected.('Side'), SIDES(ss));
                second_filter = is_freq & is_timepoint(noise_filter) & is_side;

                t_subset = t(second_filter);
                y_subset = y(second_filter);
                scatter(t_subset, y_subset, [], colors(ii), SIDES_SYMBOLS(ss))
            end
        end
        hold off
        legend(NAME_TIMEPOINTS_SIDES, 'Location', 'southeast', 'NumColumns', 2)
%         xlabel('Date')
%         ylabel('Threshold (dB)')
        ylim([10, 90])
        title([CONTROL_NOISE_LABEL{nn}, ' ', num2str(FREQ(ff)), ' Hz'])
    end

    ht.TileSpacing = 'tight';
    ht.Padding = 'tight';
    %     plot_title = ['ABR @ ', num2str(this_freq), ' Hz'];
    %     title(t, plot_title)
    %     ylabel(t, 'Amplitude (nV)')

    xlabel(ht, 'Date')
    ylabel(ht, 'DPOAE Threshold (dB)')
    % Maximize figure window size
    fig.WindowState = 'maximized';
end

%% Plot level functions
for mm = 1:n_metrics % select metric type
    for nn =  1:2 % nn = 1 is control, nn=2 is noise
        for gg = 1:N_GROUPS 
            %get group
            is_group = any(finalTable_metadata.('Batch') == BATCHGROUPS{gg}, 2);

            n_pts = zeros(n_freq, n_timepoints);
            for fff = 1:n_freq
                % one figure per frequency
                isfreq = finalTable_metadata.('Frequency')==FREQ(fff);

                figure
                hold on
                for ii = 1:n_timepoints
                    data_selected = finalTable_metadata(is_group & isfreq & metric_selector{mm} & control_noise_selector{nn} & timepoints{ii}, :);
                    x = data_selected.Properties.VariableNames(COLS_METRICVALS);
                    y = data_selected(:, COLS_METRICVALS);

                    n_pts(fff, ii) = sum(isfreq);
                    ymean = nanmean(y{:,:}, 1);
                    yste = nanstd(y{:,:}, 1)/sqrt(n_pts(fff, ii));

                    % sort points in order of ascending x value to avoid criss
                    % crossing of plot lines.
                    x_double = cellfun(@str2double, x); % convert cell array of char vectors to double matrix
                    [sortedX, sortIndex] = sort(x_double);
                    sortedYmean = ymean(sortIndex);
                    sortedYstd = yste(sortIndex);

                    % Remove nan values to plot uninterrupted curve
                    idx = ~(isnan(sortedYmean));

                    errorbar(sortedX(idx), sortedYmean(idx), sortedYstd(idx), colors(ii), 'LineWidth', LINEWIDTH)
                end
                hold off
                legend(NAME_TIMEPOINTS, 'Location','northwest')
                xlabel('Amplitude (dB)')
                ylabel(YLABEL_METRIC{mm})
                ylim(YLIM_METRIC{mm})
                xlim([-30, 100])
                title([METRIC_NAMES{mm}, ' ', CONTROL_NOISE_LABEL{nn}, ', group ', GROUP_NAMES{gg}, ', ', num2str(FREQ(fff)), ' Hz. N=', num2str(n_pts(fff, 1))]);
                set(gca,'FontSize', FONTSIZE) % Set text size of current axis

                % Save figure to results folder
                save_filename = ['levelfunction_', METRIC_NAMES{mm}, '_', CONTROL_NOISE_LABEL{nn}, '_group', GROUP_NAMES{gg}, '_freq', num2str(FREQ(fff))];
                savefig_multiformat(gcf, SAVE_PATH, save_filename)
            end
        end
    end
end