% edit_amplitude_singleABRcsv.m
%
% Read single-trace ABR response data from folder containing CSV
% spreadsheet data, and edit dB level of CSV file.
%
% To use: run file, then select folder containing CSV files in pop-up
% dialog box.
%
% Dependencies: natsortorder.m, import_ABRcsv.m
% Last edit: 8-24-21
%
% Author: George Liu

% %% Get folder path
% selpath = uigetdir();
% dir_csv = dir(fullfile(selpath, '*.csv'));
% 
% %% Import data
% A_length = length(dir_csv);
% 
% %Obtain natural order of file names
% filenames = cell(A_length, 1);
% for i = 1:A_length
%     filenames{i} = dir_csv(i).name;
% end
% filename_natorder = natsortfiles(filenames);

%%
[filename, path] = uigetfile('*.csv');
% [X_csv, A_csv, freq_csv] = import_ABRcsv(filename, path);
% disp(A_csv)

%%
% opts = spreadsheetImportOptions;
% opts = detectImportOptions(fullfile(path, filename), 'VariableNamingRule', 'preserve', 'ReadVariableNames', false);
opts = detectImportOptions(fullfile(path, filename), 'FileType', 'text', 'ReadVariableNames', false);
opts.VariableTypes{4} = 'string';
opts.VariableTypes{7} = 'string';
% opts = setvartype(opts,{'subject', 'memo'}, {'text'});
opts.VariableNamingRule = 'preserve';
this_table = readtable(fullfile(path, filename), opts);
this_table = readtable(fullfile(path, filename), 'VariableNamingRule', 'preserve', );
A_csv = this_table{9,2};

A_revised = 90-(2*(90-A_csv));

this_table{9,1} = A_revised;
outfilename = '20210814_b6m4_abr_right_try3.csv';
writetable(this_table, outfile_name)
