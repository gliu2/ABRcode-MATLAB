function [peak_pt, trough_pt, amp, latency] = get_wave1_averageABR(x, y)
% Measure wave 1 in average ABR waveform using peak detection algorithm.
% Calculates peak-peak wave 1 amplitude from wave 1 peak to following
% trough.
%
% Output: 
%
% 7/19/2021 George Liu
% Dependencies: PTDetect.m

% Constants
PEAK_SENSITIVITY = 0.1;
PEAK_START = 1; % start time (ms) for searching for wave 1 peak
PEAK_END = 2;
TROUGH_START = 1.5; % ms
TROUGH_END = 2.5; 

% get wave 1 peak coordinates
[P, T] = PTDetect(y, PEAK_SENSITIVITY);

x_peak = x(P);
y_peak = y(P);

is_ind_peak_candidates = (x_peak > PEAK_START) & (x_peak < PEAK_END);
if any(is_ind_peak_candidates)
    [peak_y, ~] = max(y_peak(is_ind_peak_candidates));
    ind = find(y_peak == peak_y);
    peak_pt = [x_peak(ind), y_peak(ind)];
else
    is_ind_peak_candidates_alt = (x > PEAK_START) & (x < PEAK_END);
    [peak_y, ~] = max(y(is_ind_peak_candidates_alt));
    ind = find(y == peak_y);
    peak_pt = [x(ind), y(ind)];
end

% get wave 1 following trough coordinates
x_trough = x(T);
y_trough = y(T);

is_ind_trough_candidates = (x_trough > TROUGH_START) & (x_trough < TROUGH_END) & (x_trough > peak_pt(1));
if ~any(is_ind_trough_candidates)
    disp('Warning: Trough not found')
    peak_pt = [NaN, NaN];
    trough_pt = [NaN, NaN];
    amp = NaN;
    latency = NaN;
    return
end
[~, ind2] = find(is_ind_trough_candidates, 1);
trough_pt = [x_trough(ind2), y_trough(ind2)];

% calculate wave 1 amplitude and latency
amp = peak_pt(2) - trough_pt(2);
latency = peak_pt(1);

end