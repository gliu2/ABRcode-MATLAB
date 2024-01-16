function y = notch_filter(x)
% Apply a 2nd order notch filter with center frequency of 1901.1 Hz, with 3 dB bandwidth of
% 200 Hz.
%
% Last edit: 8/22/23 George Liu
%
% Dependencies: none

SAMPLE_PERIOD = 5.12; % sample period of 5.12 us/sample in single trace
SAMPLE_PERIOD_MS = SAMPLE_PERIOD/1000; % sample period in ms
Fs = 1000/SAMPLE_PERIOD_MS; % Sampling frequency = 195312.50 Hz

% create notch filter with center frequency of 1901.1 Hz, with 3 dB bandwidth of
% 300 Hz.
d = fdesign.notch('N,F0,BW', 4, 1901.1, 300, Fs); 
filt = design(d,'Systemobject',true);

y = filt(x);

end