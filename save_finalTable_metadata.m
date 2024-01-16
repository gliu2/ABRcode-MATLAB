function finalTable_metadata = save_finalTable_metadata(varargin)
% Aggregate complete set of ABR single and average trace data with 
% mouse metadata from all mouse experiments.
%
% Load data from both ears (left and right). Includes single trace
% distribution of variables (eg wave 1 amp).
%
% 9/11/2022 George Liu
%
% Dependencies: import_averageABRcsv.m, savefig_multiformat.m,
% get_wave1_averageABR.m, get_roc_innerprod_ABR,
% get_thresh_averageABR_liberman.m, get_thresh_averageABR_oghalai.m

%% Constants
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving variable to .mat file
% ABR_PATH = '\\Ricci-abr\d\George\ABR'; % path to data saved on ABR computer - slow to load in data
ABR_PATH = 'D:\George-abr\ABR'; % path to local copy of data, 12-28-21

%% Load and analyze ABR data - calculate metrics and thresholds
% Initialize cache variables
close all
final_table = [];
    
%% Load data from all ABR analysis output table files into one aggregate
% table. Iterate by experimental date
fileList = dir(ABR_PATH);
alldates = {fileList(3:end).name}; % start at 3 to skip . and ..
n_dates = length(alldates);
for i=20:n_dates
    this_date = alldates{i};
    disp(['Working on experiment date ', num2str(i), ' out of ', num2str(n_dates), ': ', this_date])
    
    fileListMice = dir(fullfile(ABR_PATH, this_date));
    these_mice = {fileListMice(3:end).name};

    % Iterate over each mouse in experiment date
    num_mice = length(these_mice);
    for j=1:num_mice
        this_date_mouse = these_mice{j};
        disp(['  Working on date_mouse experiment ', num2str(j), ' out of ', num2str(num_mice), ': ', this_date_mouse])

        % Get path to ABR data for this mouse on this experimental date
        this_path = fullfile(ABR_PATH, this_date, this_date_mouse, 'analyze');
        left_right_abr_csv_files = dir(fullfile(this_path, '*.csv'));
        left_right_abr_csv_filenames = {left_right_abr_csv_files(:).name};
        num_left_right_abr_csv_files = length(left_right_abr_csv_files);
        
        % Iterate over left and right ABR data
        for k=1:num_left_right_abr_csv_files
            if strcmp(this_date_mouse, '20210807_b5m5')
                aba = 1+1;
            end
            filename = left_right_abr_csv_filenames{k};

            % Analyze ABR
            output_table = plotstack_averageABR(this_path, filename);

            % 11-4-21 added
            % Ensure number of stimulus levels (dB) are same in output table
            % and cached table results to allow vertical concatenation
            if ~isempty(final_table)
                stimulus_levels_output = output_table.Metric_values(1,:);
                stimulus_levels_cache = final_table.Metric_values(1,:);

                % check if stimulus levels are same
                num_levels_output = numel(stimulus_levels_output);
                num_levels_cache = numel(stimulus_levels_cache);
                if num_levels_output ~= num_levels_cache
                    % Merge stimulus levels of two tables
                    % cache table
                    outstanding_levels_output = setdiff(stimulus_levels_output, stimulus_levels_cache);
                    if ~isempty(outstanding_levels_output)
                        % Insert missing levels
                        add_cols = [];
                        for z=outstanding_levels_output
                            new_col = NaN(size(final_table.Metric_values, 1), 1);
                            new_col(1:3) = z;
                            add_cols = [add_cols, new_col];
                        end
                        new_final_table_metricvals = [final_table.Metric_values, add_cols];
                        [new_stimulus_levels_cache, ind_levels_sorted] = sort(new_final_table_metricvals(1,:), 'descend');
                        new_final_table_metricvals = new_final_table_metricvals(:, ind_levels_sorted);
                        final_table.Metric_values = new_final_table_metricvals;
                    end

                    % output table
                    outstanding_levels_cache = setdiff(stimulus_levels_cache, stimulus_levels_output);
                    if ~isempty(outstanding_levels_cache)
                        % Insert missing levels
                        add_cols = [];
                        for z=outstanding_levels_cache
                            new_col = NaN(size(output_table.Metric_values, 1), 1);
                            new_col(1:3) = z;
                            add_cols = [add_cols, new_col];
                        end
                        new_output_table_metricvals = [output_table.Metric_values, add_cols];
                        [new_stimulus_levels_output, ind_levels_sorted] = sort(new_output_table_metricvals(1,:), 'descend');
                        new_output_table_metricvals = new_output_table_metricvals(:, ind_levels_sorted);
                        output_table.Metric_values = new_output_table_metricvals;
                    end
                end
            end

            % Cache table results
            final_table = vertcat(final_table, output_table);
        end
    end
    
    % Checkpoint final_table variable in case runtime gets interrupted
    save(fullfile(SAVE_PATH, 'final_table.mat'), 'final_table')
end

%% Add metadata to table
finalTable = final_table;
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

end