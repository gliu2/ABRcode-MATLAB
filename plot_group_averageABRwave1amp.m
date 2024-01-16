% Wrapper script to plot average wave 1 amplitude pre, post24hour, and post1week 
% Group mice by noise exposure group
% Code adapted from plot_ABRthreshold_levelfunctions_onemouse.m
%
% SAVE_PATH automatically set to "Results" folder to save output files
%
% Dependencies: plotstack_averageABR.m
%
% 12/8/21 George Liu

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
YLIM_METRIC = {[0, 6250]};
YLABEL_METRIC = {'Wave 1 amplitude (nV)'}; 

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
            if strcmp(SIDES{nn}, "both")
                is_side = ones(size(finalTable_metadata.Side));
            else
                is_side = finalTable_metadata.Side == SIDES{nn};
            end
            
            fig = figure;
            t = tiledlayout(n_conditions, n_freq, 'TileIndexing', 'rowmajor');
            
            for cc = 1:n_conditions % Iterate over control then noise mice
                n_pts = zeros(n_freq, n_timepoints);
                for fff = 1:n_freq
                    % one figure per frequency
                    isfreq = finalTable_metadata.('Frequency')==FREQ(fff);

                    nexttile(t)
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

        %                 % 9-6-22: try normalizing wave1amp level function of single
        %                 % mouse by max stimulus amplitude to see if control mouse then
        %                 % has no change in wave 1 amp at 1 week
        %                 sortedYmean_real = sortedYmean_real/max(sortedYmean_real);

                        % Remove stimulus levels that are not divisible by 10
                        isdiv10 = mod(sortedX_real, 10)==0; 
                        errorbar(sortedX_real(isdiv10), sortedYmean_real(isdiv10), sortedYstd_real(isdiv10), colors(ii), 'LineWidth', LINEWIDTH)
        %                 plot(sortedX_real, sortedYmean_real, colors(ii), 'LineWidth', LINEWIDTH)

                        % Determine if pre and post are significantly different
                        % at each stimulus level
                        if ii==1
                            y_pre = y{:,:};
                            y_pre = y_pre(:, sortIndex);
                        else
                            y_post = y{:,:};
                            y_post = y_post(:, sortIndex);
                            for zz = 1:numel(sortedX) % iterate thru stimulus levels
                                this_prey = y_pre(:, zz);
                                this_posty = y_post(:, zz);

                                if all(isnan(this_prey)) || all(isnan(this_posty))
                                    continue
                                end

                                % Remove NaN values
                                this_prey = this_prey(~isnan(this_prey));
                                this_posty = this_posty(~isnan(this_posty));

                                % Compare with t-test
                                is_different = ttest2(this_prey, this_posty, 'Tail', 'both');
                                if is_different && mod(sortedX(zz), 10)==0
                                    % asterisk on wave 1 amp curve if significantly different from baseline
                                    text(sortedX(zz), mean(this_posty), '*', ...
                                        'Color', colors(ii), 'FontSize', FONTSIZE*2, 'FontWeight', 'bold') 
                                end
                            end
                        end

                    end

                    %TODO: add statistical significance test between pre and post
                    %noise curves

                    hold off
                    legend(NAME_TIMEPOINTS, 'Location','northwest')

                    xlabel('Amplitude (dB)')
                    ylabel(YLABEL_METRIC{mm})
                    ylim(YLIM_METRIC{mm})
        %             ylim([-0.1, 1.1]) % For normalized wave1amp
                    xlim([-30, 100])
            %         title([METRIC_NAMES{mm}, ' ', ', ', num2str(FREQ(fff)), ' Hz. N=', num2str(n_pts(fff, 1))]); % TODO: fix n_pts, which is currently 33 times (n_timepts*num_levels) the actual number of mice
                    title([CONTROL_NOISE_LABEL{cc}, ', ', num2str(FREQ(fff)), ' Hz. N=', num2str(n_pts(fff, ii)), ' ears'], 'Interpreter','none');
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

disp('Done.')