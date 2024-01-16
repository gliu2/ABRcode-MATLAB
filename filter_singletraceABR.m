function X_csv_new = filter_singletraceABR(X_csv, N, varargin)
% Filters single traces in single trace ABR file. Uses Blackman filter.
%
% Stopped using because this implementation, with MATLAB's filter function,
% shifts the curve by an amount proportional to the filter order.
%
% Input: X_csv - ABR single trace dataset.
%            Matrix of size (SAMPLES, m_traces), for set of single ABR
%            traces at same dB level. Columns are individual trace data.
%
%        N - order of filter
%
%        hp_f - high pass frequency cutoff (Hz)
%
%        lp_f - low pass frequency cutoff (Hz)
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

if nargin == 3
    hp_f = varargin{1};
elseif nargin == 4
    hp_f = varargin{1};
    lp_f = varargin{2};
else
    error('Error, wrong number of inputs to filter_singletraceABR.m')
end

Fs = 195300; % Sampling rate: 195.3 kHz

% Bandpass Blackman frequency
Fn = Fs/2;  % Nyquist Frequency
N = 2*fix(N/2); % N must be even
wndw = blackman(N+1);
if nargin == 3
    b = fir1(N, hp_f/Fn, 'high', wndw);
elseif nargin == 4
    b = fir1(N, [hp_f/Fn, lp_f/Fn], 'bandpass', wndw);
end

 
X_csv_new = filter(b, 1, X_csv);

end