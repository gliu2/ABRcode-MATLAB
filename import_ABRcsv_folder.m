% import_ABRcsv_folder.m
%
% Read single-trace ABR response data from folder containing CSV
% spreadsheet data, one CSV per dB level.
%
% To use: run file, then select folder containing CSV files in pop-up
% dialog box.
%
% Dependencies: natsortorder.m
% Last edit: 6-8-19
%
% Author: George Liu

%% Get folder path
selpath = uigetdir();
dir_csv = dir(fullfile(selpath, '*.csv'));

%% Import data
A_length = length(dir_csv);
X_csv = cell(A_length, 1);
A_csv = zeros(A_length, 1);
freq_csv = zeros(A_length, 1);

%Obtain natural order of file names
filenames = cell(A_length, 1);
for i = 1:A_length
    filenames{i} = dir_csv(i).name;
end
filename_natorder = natsortfiles(filenames);

for i = 1:A_length
%     filename = dir_csv(i).name;
%     [X_csv{i}, A_csv(i), freq_csv(i)] = import_ABRcsv(filename, selpath);
    [X_csv{i}, A_csv(i), freq_csv(i)] = import_ABRcsv(filename_natorder{i}, selpath);
end
