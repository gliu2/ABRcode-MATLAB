%% import_ABRcsv.m
%
% Read ABR response for CSV spreadsheet
% Data values are: x-axis time (ms); y-axis voltage (nV); sampling rate 200K Hz
% 2/28/2019 George Liu
% Dependencies: none

function [X_csv, A_csv, freq_csv] = import_ABRcsv(filename, path)
% [filename,path] = uigetfile('*.csv');
filepath = fullfile(path, filename);

% Row column offsets for cell A3 (top left corner)
R1 = 2;
C1 = 0;
M = csvread(filepath,R1,C1); % m x SAMPLES matrix

% Remove last column if all zeros
% When the csvread function reads data files with lines that end with a nonspace delimiter, such as a semicolon, it returns a matrix, M, that has an additional last column of zeros.
if any(M(:,end))==0
    M = M(:,1:end-1);
end

% Get ABR input frequency and amplitude
R2 = 1;
C2 = 7; % H2 Frequency
C3 = 8; % I2 input dB amplitude
M2 = csvread(filepath,R2,C2, [R2 C2 R2 C3]);
freq_csv = M2(1); % Hz
A_csv = M2(2); % dB

X_csv = M'; % SAMPLES x m matrix
% M_avg = mean(X_csv, 2); % SAMPLES x 1 matrix

end