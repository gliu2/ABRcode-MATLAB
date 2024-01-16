% load_all_ABR_DPOAE_analysis
%
% Load all ABR or DPOAE analysis files in results folder.
%
% Adapted from plot_ABRthreshold_all_scatter.m
%
% 11/11/21 George Liu
% Dependencies: none

function finalTable = load_all_ABR_DPOAE_analysis(load_path, is_abr)

if is_abr
    fileList = dir(fullfile(load_path, '*_abr_output.xlsx'));
else
    fileList = dir(fullfile(load_path, '*_dpoae_output.xlsx'));
end
n_files = length(fileList);
for j=1:n_files
    disp(['Working on file ', num2str(j), ' out of ', num2str(n_files), ': ', fileList(j).name])
    outputFile = fullfile(load_path, fileList(j).name);
    outputTable = readtable(outputFile);
    
    if is_abr
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
    end
    
    % Concatenate tables
    if j==1
        finalTable = outputTable;
    else
        finalTable = outerjoin(finalTable,outputTable,'MergeKeys', true);
    end
end