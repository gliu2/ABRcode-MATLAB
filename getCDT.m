function [f1, f1_amp, f2, f2_amp, cdt_freq, cdt_amp, z_score] = getCDT(A, varargin)
%GETCDT - Calculate cubic distortion tone from DPOAE measurement
%   Returns calculated fundamental frequenceis and their 2f1-f2 cdt. 
%   Assumes f1 and f2 generate largest two peaks in DPOAE.
%   Based on https://github.com/CDTbot/CDTbot
%
%   Optional inputs: getCDT(A, f1, f2)
% 
%   7/15/21
%   George Liu

FREQ_INTERVAL = 97656.25/2048; % Hz

% Fundamental freqencies
[~, I] = sort(A, 'descend');
if nargin == 1
    f1 = (I(1) - 1) * FREQ_INTERVAL;
    f2 = (I(2) - 1) * FREQ_INTERVAL;
    
elseif nargin == 3
    f1 = varargin{1};
    f2 = varargin{2};
else
    disp('Warning: Wrong number of inputs to getCDT!')
    f1 = (I(1) - 1) * FREQ_INTERVAL;
    f2 = (I(2) - 1) * FREQ_INTERVAL;
end

% Ensure f2 > f1
if f1 > f2
    temp = f1;
    f1 = f2;
    f2 = temp;
end

% Amplitudes of fundamental frequencies
f1_ind = round(f1/FREQ_INTERVAL) + 1;
f2_ind = round(f2/FREQ_INTERVAL) + 1;
f1_amp = A(f1_ind);
f2_amp = A(f2_ind);

% Calculate cubic distortion tone 
cdt_freq = 2*f1 - f2;
if cdt_freq < 0
    cdt_freq = 0;
end
x = FREQ_INTERVAL * (0:(length(A) - 1));
[~, index_cdt] = min(abs(x - cdt_freq));
cdt_amp = A(index_cdt);

% Calculate noise floor average and standard deviation.
% Uses 10 values from either side of the CDT, not including the 2 values 
% directly adjacent to the DP.
noise_ind_min = index_cdt - 11;
noise_ind_max = index_cdt + 11;
if noise_ind_min >0
    noise_ind = [noise_ind_min:(index_cdt - 2), (index_cdt+2):noise_ind_max];
else
    noise_ind = length(A) - 100:length(A)-10;
end
noise_ave = mean(A(noise_ind));
noise_std = std(A(noise_ind));
z_score = (cdt_amp - noise_ave)/noise_std;
end

