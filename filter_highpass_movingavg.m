function X_csv_new = filter_highpass_movingavg(X_csv, N)
% Filters single traces in single trace ABR file. Uses moverage average filter.
%
% Input: X_csv - ABR single trace dataset.
%            Matrix of size (SAMPLES, m_traces), for set of single ABR
%            traces at same dB level. Columns are individual trace data.
%
%        N - order of filter
%
% Output: X_csv_new - ABR single trace dataset, after filtering.
%            Matrix of size (SAMPLES, m_traces/2), for set of single ABR
%            traces at same dB level. Columns are individual trace data.
%
%
% Dependencies: none
% George Liu 
%
% Last edit 7/19/23

Fs = 195300; % Sampling rate: 195.3 kHz

% Moving average window
windowWidth = N; % Whatever you want.
kernel = ones(windowWidth,1) / windowWidth;

% This implementation of filtering, rather than filter function, doesn't
% translate filtered signal.
m_traces = size(X_csv, 2);
out = zeros(size(X_csv));
for i = 1:m_traces
    out(:, i) = conv(X_csv(:, i), kernel, 'same');
end
% X_csv_new = out; % low pass output
X_csv_new = X_csv - out; % high pass output

% lowPassed = filter(kernel, 1, X_csv);
% X_csv_new = X_csv - lowPassed; % high pass output

end