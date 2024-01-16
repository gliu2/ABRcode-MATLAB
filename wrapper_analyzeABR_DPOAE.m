% Wrapper script to run analysis on all ABR or DPOAE csv files in folder
%
% **BEFORE RUNNINNG**: MANUALLY SET IS_ABR CONSTANT TO 0 (ABR) OR 1 (DPOAE)TO SELECT FILE
% TYPE TO ANALYZE
%
% SAVE_PATH automatically set to "Results" folder to save output files
%
% 7/15/21 George Liu

close all

%% Constants
IS_ABR = 1; % MANUALLY SET TO 0 (ABR) OR 1 (DPOAE)
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures

% Initialize cache variables
final_table = [];

%% Select folder with all subfolders containing CSV files
path = uigetdir();
files = dir(path);
dirFlags = [files.isdir] & ~strcmp({files.name},'.') & ~strcmp({files.name},'..');
subfolders = files(dirFlags);

for j=1:length(subfolders)
    disp(['Working on folder ', num2str(j), ' out of ', num2str(length(subfolders)), ': ', subfolders(j).name])
    if IS_ABR
        this_path = fullfile(path, subfolders(j).name, 'analyze');
    else
        this_path = fullfile(path, subfolders(j).name);
    end
    fileList = dir(fullfile(this_path, '*.csv'));
    n_files = length(fileList);

    %% Analyze ABR
    for i=1:n_files
        filename = fileList(i).name;
        
        if IS_ABR
            if contains(filename, "mask", 'IgnoreCase', true)
                % Masked ABR
                % TODO: add masked ABR analysis function
                continue
            else
                % ABR (unmasked)
                output_table = plotstack_averageABR(this_path, filename);
            end
        else
            % DPOAE
            output_table = plotstack_DPOAE(this_path, filename);
        end
        
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

% Write combine output results to excel sheet
% [~, baseFileNameNoExt, ~] = fileparts(filename);
k = strfind(path, '2021');
kk = strfind(path, filesep);
is_postslash = kk>k;
if any(is_postslash)
    ind_postslash = kk(is_postslash);
    baseFileNameNoExt = path(k:ind_postslash(1)-1);
else
    baseFileNameNoExt = path(k:end);
end

if IS_ABR
    file_ext = '_abr_output.xlsx';
else
    file_ext = '_dpoae_output.xlsx';
end
table_filename = [baseFileNameNoExt, file_ext];
writetable(final_table, fullfile(SAVE_PATH, table_filename), 'Sheet', 1, 'Range', 'A1')

disp('Done.')