%% import_ABRcsv.m
%
% Read single trace ABR response for CSV spreadsheet
% Data values in CSV file are: x-axis time (samples); y-axis voltage (V); sampling rate 200K Hz
%
% Sampling rate: 195.3 kHz
% Sampling period: 5.12 us/sample
% 1952 samples (8x more sampling frequency than average CSV trace)
% Total time: 9.9942 ms (first sample is at 1 sample period, not at 0 ms)
%
% % Output: 
%       X_csv    - SAMPLES x m matrix of m single traces, each with SAMPLES number of data samples. In
%                this case there are 1952 data samples. Voltage units are converted to nV (originally in V).
%
%       freq_csv - frequency of this set of single. Scalar
%       A_csv    - stimulus level of this set of single traces. Scalar
%
% George Liu
% Dependencies: none
%
% Last edit: 
% 12/29/21 - replaced csvread with readmatrix to load data - slightly faster. 
% Kept second call to csvread to read in stimulus amplitude (dB) and
% frequency.
%
% 12/10/21 
% - edit output to be in units  of nV (original units are V)
% - remove DC shift from single traces, since it is not removed by TDT
% software (but does remove for average trace)
% Created 2/28/2019 

function [X_csv, A_csv, freq_csv] = import_ABRcsv(filename, path)
% Initialize constants
V_TO_NV_CONVERSIONFACTOR = 10^9; % multiplication factor to convert original units (V) to nV

% [filename,path] = uigetfile('*.csv');
filepath = fullfile(path, filename);

% Row column offsets for cell A3 (top left corner)
% R1 = 2;
% C1 = 0;
% M = csvread(filepath,R1,C1); % m x SAMPLES matrix
% 
% % Remove last column if all zeros
% % When the csvread function reads data files with lines that end with a nonspace delimiter, such as a semicolon, it returns a matrix, M, that has an additional last column of zeros.
% if any(M(:,end))==0
%     M = M(:,1:end-1);
% end
M = readmatrix(filepath, 'Range', [3,1]);

% Get ABR input frequency and amplitude
R2 = 1;
C2 = 7; % H2 Frequency
C3 = 8; % I2 input dB amplitude
M2 = csvread(filepath,R2,C2, [R2 C2 R2 C3]);

% M2 = readmatrix(filepath, 'Range', 
freq_csv = M2(1); % Hz
A_csv = M2(2); % dB

X_csv = M'; % SAMPLES x m matrix
% M_avg = mean(X_csv, 2); % SAMPLES x 1 matrix

% 12/10/21 - convert units V to nV
X_csv = X_csv*V_TO_NV_CONVERSIONFACTOR; 
dc_shift = mean(X_csv, 1);
X_csv = X_csv - dc_shift; % Remove DC shift from single traces

end