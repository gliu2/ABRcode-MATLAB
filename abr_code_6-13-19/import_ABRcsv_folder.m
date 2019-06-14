% import_ABRcsv_folder.m
%
% Read single-trace ABR response data from folder containing CSV
% spreadsheet data, one CSV per dB level.
%
% To use: run file, then select folder containing CSV files in pop-up
% dialog box.
%
% IMPORTANT: To use full time of ABR traces, comment out last cell.
% To remove first 2 ms of each ABR single trace, uncomment last cell.
%
% Dependencies: natsortorder.m, import_ABRcsv.m
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
    [X_csv{i}, A_csv(i), freq_csv(i)] = import_ABRcsv(filename_natorder{i}, selpath);
end

%% 6-12-19: Obtain only a fixed number of single traces per dB level
IMPORT_NUM_TRACES = 514;
for i = 1:A_length
    X_csv{i} = X_csv{i}(:, 1:IMPORT_NUM_TRACES);
end

%% Truncate first 2 ms of (noise) data if tone pip starts at t=2 ms (manually comment out for traces with tone pip at 0 ms)
SAMPLES = size(X_csv{1}, 1);
SAMPLING_RATE = 200000; % Hz
dt = 1/SAMPLING_RATE * 1000; % ms
T_total = (SAMPLES-1)*dt; % ms
x = 0:dt:T_total;

TIME_TONEPIP = 2; % ms
ind_tonepip = find(x==TIME_TONEPIP); 
for i = 1:A_length
    X_csv{i} = X_csv{i}(ind_tonepip:end, :);
end

%% 6-12-19: Ensure trace is 7.755 ms long, even if tone pip is played at start
SAMPLES = size(X_csv{1}, 1);
SAMPLING_RATE = 200000; % Hz
dt = 1/SAMPLING_RATE * 1000; % ms
T_total = (SAMPLES-1)*dt; % ms
x = 0:dt:T_total;

NEW_T_TOTAL = 7.755; % ms
tol = 1e-3; 
ind_Tend = find(abs(x-NEW_T_TOTAL)<tol);
for i = 1:A_length
    X_csv{i} = X_csv{i}(1:ind_Tend, :);
end
