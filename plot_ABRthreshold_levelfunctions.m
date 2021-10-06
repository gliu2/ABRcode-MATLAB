% plot_ABRthreshold_levelfunctions.m
%
% 9/2/21 George Liu
% Dependencies: get_mousefile_metadata.m, savefig_multiformat.m

%% Constants
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures
LOAD_PATH = 'd:\users\admin\Documents\George\Results_9-1-21';
BATCHGROUPS = {1:4, 5:6, 7};
N_GROUPS = numel(BATCHGROUPS);
METRIC_NAMES = {'Liberman', 'Oghalai', 'Innerprod', 'Innerprod-AUC'};
n_metrics = numel(METRIC_NAMES);
FREQ = [8000, 16000, 32000];
n_freq = numel(FREQ);
COLS_METRICVALS = 22:size(finalTable_metadata, 2); % column numbers of metric values columns in finalTable_metadata table
CONTROL_NOISE_LABEL = {'control', 'noise'};
YLIM_METRIC = {[-0.4, 1], [0, 6250], [0, 2.5*10^-9], [0.5, 1]};
YLABEL_METRIC = {'Correlation coefficient (a.u.)', 'Max p-p amplitude (nV)', 'Mean inner product (nv^2)', 'AUC (a.u.)'}; 

% Initialize groups of measurements
% Frequency = freq_unique;
% Stimulus_amplitude = zeros(n_freq, num_levels);
% Threshold_liberman = zeros(n_freq, 1);
% Threshold_oghalai = zeros(n_freq, 1);
% Threshold_innerprod = zeros(n_freq, 1);
% Metric_liberman = zeros(n_freq, num_levels);
% Metric_oghalai = zeros(n_freq, num_levels);
% Metric_innerprod = zeros(n_freq, num_levels);
% Metric_innerprod_auc = zeros(n_freq, num_levels);

%% Load data from all ABR analysis output table files into one aggregate table
fileList = dir(fullfile(LOAD_PATH, '*_abr_output.xlsx'));
n_files = length(fileList);
for j=1:n_files
    disp(['Working on file ', num2str(j), ' out of ', num2str(n_files), ': ', fileList(j).name])
    outputFile = fullfile(LOAD_PATH, fileList(j).name);
    outputTable = readtable(outputFile);
    
    % Rename Metric_values_n variable names to stimulus amplitudes
    variablenames = outputTable.Properties.VariableNames;
    n_variables = numel(variablenames);
    new_variablenames = variablenames;
    for i=1:n_variables
        if contains(variablenames{i}, 'Metric_values', 'IgnoreCase',true)
            stimulus_amplitude = outputTable.(variablenames{i})(1); % stimulus amplitude
            new_variablenames{i} = num2str(stimulus_amplitude);
        end
    end
    outputTable.Properties.VariableNames = new_variablenames;
    
    % Concatenate tables
    if j==1
        finalTable = outputTable;
    else
        finalTable = outerjoin(finalTable,outputTable,'MergeKeys', true);
    end
end

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
% Metric filters
is_liberman = contains(finalTable_metadata.('Metric'), 'Liberman', 'IgnoreCase',true);
is_oghalai = contains(finalTable_metadata.('Metric'), 'Oghalai', 'IgnoreCase',true);
is_innerprod_auc = contains(finalTable_metadata.('Metric'), 'Innerprod_auc', 'IgnoreCase',true);
is_innerprod = contains(finalTable_metadata.('Metric'), 'Innerprod', 'IgnoreCase',true) & ~is_innerprod_auc;
metric_selector = {is_liberman, is_oghalai, is_innerprod, is_innerprod_auc};

% Noise control filter
is_noise = finalTable_metadata.('Is_noise_exposed');
is_control = ~is_noise;
control_noise_selector = {is_control, is_noise};

% Timepoint
is_pre = strcmp(finalTable_metadata.('Date_1'), finalTable_metadata.('Date'));
is_post24 = strcmp(finalTable_metadata.('Date_2'), finalTable_metadata.('Date'));
is_post1w = strcmp(finalTable_metadata.('Date_3'), finalTable_metadata.('Date'));
timepoints = {is_pre, is_post24, is_post1w};
n_timepoints = numel(timepoints);
name_timepoints = ["Pre", "Post 24h", "Post 1w"];

%% Plot threshold shift per frequency
colors = ['b', 'g', 'r'];

for gg = 1:N_GROUPS
    %get group
    is_group = any(finalTable_metadata.('Batch') == BATCHGROUPS{gg}, 2);
    
    figure
    hold on
    for ii = 1:n_timepoints
        data_selected = finalTable_metadata(is_group & is_liberman & is_noise & timepoints{ii}, :);
        x = data_selected.('Frequency');
        y = data_selected.('Threshold');

        ymean = zeros(n_freq, 1);
        ystd = zeros(n_freq, 1);
        n_pts = zeros(n_freq, 1);
        for ff = 1:n_freq
            n_pts(ff) = sum(x==FREQ(ff));
            ymean(ff) = nanmean(y(x==FREQ(ff)));
            ystd(ff) = nanstd(y(x==FREQ(ff)))/sqrt(n_pts(ff));
        end

    %     scatter(x, y)
        errorbar(FREQ, ymean, ystd, colors(ii))
    end
    hold off
    legend(name_timepoints, 'Location','northwest')
    xlabel('Frequency (Hz)')
    ylabel('Threshold (dB)')
    ylim([-10, 90])
    title([METRIC_NAMES{1}, ' noise, group ', num2str(gg)]);
    disp(n_pts)
end

%% Plot level functions
for mm = 1:n_metrics % select metric type
    for nn =  1:2 % nn = 1 is control, nn=2 is noise
        for gg = 1:N_GROUPS 
            %get group
            is_group = any(finalTable_metadata.('Batch') == BATCHGROUPS{gg}, 2);

            for fff = 1:n_freq
                % one figure per frequency
                isfreq = finalTable_metadata.('Frequency')==FREQ(fff);

                figure
                hold on
                for ii = 1:n_timepoints
                    data_selected = finalTable_metadata(is_group & isfreq & metric_selector{mm} & control_noise_selector{nn} & timepoints{ii}, :);
                    x = data_selected.Properties.VariableNames(COLS_METRICVALS);
                    y = data_selected(:, COLS_METRICVALS);

                    n_pts = sum(isfreq);
                    ymean = nanmean(y{:,:}, 1);
                    ystd = nanstd(y{:,:}, 1)/sqrt(n_pts);

                    % sort points in order of ascending x value to avoid criss
                    % crossing of plot lines.
                    x_double = cellfun(@str2double, x); % convert cell array of char vectors to double matrix
                    [sortedX, sortIndex] = sort(x_double);
                    sortedYmean = ymean(sortIndex);
                    sortedYstd = ystd(sortIndex);

                    % Remove nan values to plot uninterrupted curve
                    idx = ~(isnan(sortedYmean));

                    errorbar(sortedX(idx), sortedYmean(idx), sortedYstd(idx), colors(ii))
                end
                hold off
                legend(name_timepoints, 'Location','northwest')
                xlabel('Amplitude (dB)')
%                 ylabel('Metric (a.u.)')
                ylabel(YLABEL_METRIC{mm})
                ylim(YLIM_METRIC{mm})
                xlim([-30, 100])
                title([METRIC_NAMES{mm}, ' ', CONTROL_NOISE_LABEL{nn}, ', group ', num2str(gg), ', ', num2str(FREQ(fff)), ' Hz']);
                disp(n_pts)

                % Save figure to results folder
                save_filename = ['levelfunction_', METRIC_NAMES{mm}, '_', CONTROL_NOISE_LABEL{nn}, '_group', num2str(gg), '_freq', num2str(FREQ(fff))];
                savefig_multiformat(gcf, SAVE_PATH, save_filename)
            end
        end
    end
end