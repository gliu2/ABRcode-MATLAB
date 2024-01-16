function [X_csv, x_ms] = import_singleABR_frequency(input_pathfile, freq, amp)
%% import_singleABR_frequency.m
%
% Read single trace ABR response at a specified frequency (and intensity) 
% using the input path of the average ABR CSV spreadsheet.
% Data values in CSV file are: x-axis time (samples); y-axis voltage (V); sampling rate 200K Hz
%
% Sampling rate: 195.3 kHz
% Sampling period: 5.12 us/sample
% 1952 samples (8x more sampling frequency than average CSV trace)
% Total time: 9.9942 ms (first sample is at 1 sample period, not at 0 ms)
%
% Inputs (optional):
%       - input_pathfile: full file path to average ABR CSV, related to single
%       ABR CSV files.
%       freq_csv - frequency of this set of single. Scalar
%       A_csv    - stimulus level of this set of single traces. Scalar
%
% % Output: 
%       X_csv    - SAMPLES x m/2 matrix of m/2 single traces, with merged 
%                  single trace pairs with alternating polarity to cancel cochlear
%                microphonic heterogeneity. Each single trace has SAMPLES number of data samples. In
%                this case there are 1952 data samples. Voltage units are
%                converted to nV (originally in V). 
%
% George Liu
% Dependencies: import_ABRcsv.m
%
% Last edit: 
% 9/21/2022

SAMPLING_RATE = 195.3; % samples / msec. Sampling rate of single trial ABR CSV file. Of note, average CSV files have different sampling rate b/c they compress number of samples 1952 -> 244.
    
[selpath_analyze, file_pre, file_ext] = fileparts(input_pathfile);
file_pattern = [file_pre, '-*', file_ext]; % include hyphen, which is more specific to single trace ABR CSV files.
selpath = fileparts(selpath_analyze); % get parent folder containing "analyze" folder

% get listing of all single trace CSV files with the file name prefix
parent_folder = fullfile(selpath, file_pattern);
dir_files = dir(parent_folder);
n_files = numel(dir_files);

% Throw warning message and abort if parent directory contains no single trace ABR CSV Files.
if isempty(dir_files) 
    disp(['Warning! No single trace ABR CSV files in parent folder: ', parent_folder])
    disp('-> Aborting single trace analysis of ABR in get_roc_innerprod_ABR.m, line 98.')
    disp('')
    return
end

for i = 1:n_files
    
    % Load next single trial ABR data (CSV file)
    this_filename = dir_files(n_files - i + 1).name;
    disp(['Working on file ', num2str(i), ' out of ', num2str(n_files), ': ', this_filename])
    [X_csv, A_csv, freq_csv] = import_ABRcsv(this_filename, selpath);
    
    if A_csv==amp && freq_csv==freq
        % merge single trace pairs with alternating polarity to cancel cochlear
        % microphonic heterogeneity
        X_csv = merge_singletraceABR_polarities(X_csv); % 1026 -> 513 single traces. Matrix of size (SAMPLES, m_traces/2)
        
        % Calculate time points of single trial ABR waveform
        num_samples = size(X_csv, 1);
        x_ms = (1/SAMPLING_RATE) * (1:num_samples); % time in ms
        
        return 
    end
end

% Return NaN if no single trace ABR is found
X_csv = NaN;
x_ms = NaN;
disp(['Warning! No single trace ABR CSV file found for input frequency ', num2str(freq), ' Hz and stimulus intensity ', num2str(amp), ' nV'])
disp('-> Returning NaN for single trace analysis of ABR in import_singleABR_frequency.m, line 75.')
