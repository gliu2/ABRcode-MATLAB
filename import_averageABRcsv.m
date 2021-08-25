function [M, A_csv, freq_csv] = import_averageABRcsv(filename, path)
%import_averageABRcsv.m - Import CSV file containing average ABR trace data
% Read average ABR response for CSV spreadsheet that is exported from
% BioSigRZ workspace.
% Data values are: x-axis time (ms); y-axis voltage (nV); sampling rate 200K Hz
%
% Output: M - n x m array of n average traces, each with m data samples. In
% this case there are 244 data samples.
%
% 7/9/2021 George Liu
% Dependencies: none

% [filename,path] = uigetfile('*.csv');
filepath = fullfile(path, filename);

% Row column offsets for cell (2, AW) 
R1 = 1; % Start at row 2
C1 = 49 - 1; % Column AW is column 49
M = csvread(filepath,R1,C1); % m x SAMPLES matrix

% Remove last column if all zeros
% When the csvread function reads data files with lines that end with a nonspace delimiter, such as a semicolon, it returns a matrix, M, that has an additional last column of zeros.
if any(M(:,end))==0
    M = M(:,1:end-1);
end

% Get ABR input frequency and amplitude
N = size(M, 1); % number of traces
R2 = 2;
R3 = R2 + N - 1;
C2 = 13; % column M index is 13, with Frequency
C3 = 14; % column N input dB amplitude
M2 = readmatrix(filepath, 'Range', [R2 C2 R3 C3]);
freq_csv = M2(:, 1); % Hz
A_csv = M2(:, 2); % dB

% Correct if headers for freq and dB are switched
headers = readmatrix(filepath, 'Range', [1 C2 1 C3], 'OutputType', 'string');
if contains(headers(1), 'level', 'IgnoreCase', true)
    freq_csv = M2(:, 2); % Hz
    A_csv = M2(:, 1); % dB
end

end