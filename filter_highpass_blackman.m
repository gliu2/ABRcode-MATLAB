function X_csv_new = filter_highpass_blackman(X_csv, windowWidth)
% Filters single traces in single trace ABR file. Uses Blackman filter.
% This implementation doesn't shift the input trace, unlike the MATLAB
% filter function.
%
% Input: X_csv - ABR single trace dataset.
%            Matrix of size (SAMPLES, m_traces), for set of single ABR
%            traces at same dB level. Columns are individual trace data.
%
%        windowWidth - order of filter
%
% Output: X_csv_new - ABR single trace dataset, after filtering.
%            Matrix of size (SAMPLES, m_traces/2), for set of single ABR
%            traces at same dB level. Columns are individual trace data.
%
%
% Dependencies: none
% George Liu 
%
% Last edit 7/28/23

% Fs = 195300; % Sampling rate: 195.3 kHz

kernel = blackman(windowWidth);
kernel = kernel / sum(kernel);

m_traces = size(X_csv, 2);
out = zeros(size(X_csv));
for i = 1:m_traces
    out(:, i) = conv(X_csv(:, i), kernel, 'same');
end
X_csv_new = X_csv - out; % High pass output
% X_csv_new = out; % Low pass output

% X_csv_new = filter(kernel, 1, X_csv);
end