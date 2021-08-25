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
                continue
            else
                % ABR (unmasked)
                output_table = plotstack_averageABR(this_path, filename);
            end
        else
            % DPOAE
            output_table = plotstack_DPOAE(this_path, filename);
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