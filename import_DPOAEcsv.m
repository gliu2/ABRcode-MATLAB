function [M, A_csv, freq_csv, f1, f2, sample_period] = import_DPOAEcsv(filename, path)
% Import CSV file containing DPOAE data
% Data values are: x-axis time (ms); y-axis (dBnV); sampling rate 200K Hz
%
% Output: M - n x m array of n average traces, each with m data samples. In
% this case there are 2047 data samples.
%
% 7/12/2021 George Liu
% Dependencies: none

% [filename,path] = uigetfile('*.csv');
filepath = fullfile(path, filename);

% Row column offsets for cell (2, AW) 
R1 = 1; % Start at row 2
C1 = 49 - 1; % Column AW is column 49
M = csvread(filepath,R1,C1); % m x SAMPLES matrix

% Remove last column if all zeros
% When the csvread function reads data files with lines that end with a nonspace delimiter, such as a semicolon, it returns a matrix, M, that has an additional last column of zeros.
while any(M(:,end))==0
    M = M(:,1:end-1);
end

% Get input frequency and amplitude
N = size(M, 1); % number of traces
R2 = 1;
R3 = R2 + N - 1;
C2 = 12; % column M index is 13, dB amplitude
C3 = 13; % column N input frequency (Hz)
C4 = 16; % column Q F1 (Hz)
C5 = 17; % column R F2 (Hz)
C6 = 44; % column AS (us/sample)
M2 = csvread(filepath,R2,C2, [R2 C2 R3 C5]);
freq_csv = M2(:, 2); % Hz
A_csv = M2(:, 1); % dB
f1 = M2(:, 5);
f2 = M2(:, 6);
sample_period = csvread(filepath,R2,C6, [R2 C6 R2 C6]);

end