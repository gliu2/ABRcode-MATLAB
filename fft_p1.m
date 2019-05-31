%% fft_p1.m
% 
% Return single-sided spectrum P1 of fast fourier transform of signal X (analyzed column by column),
% and its frequency domain f. Maxixum frequency is one half the sampling
% rate Fs.
%
% Input: X - Signal data. Matrix of column vectors. Size (rows, columns)
%        Fs - sampling rate (Hz)
%
% Output: f - frequency domain (Hz) for plotting on x-axis. Row vector.
%         P1 - single-sided power spectrum of fft. Matrix of size (rows/2+1, columns)
%
% Dependencies: none
% Last edit: 5/30/2019
%
% Author: George Liu

function [f, P1] = fft_p1(X, Fs)

% Compute the Fourier transform of the signal.
Y = fft(X); % take fft of column vectors in matrix X

% Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
L = size(X, 1);
P2 = abs(Y/L);
P1 = P2(1:L/2+1, :);
P1(2:end-1, :) = 2*P1(2:end-1, :);

% Define the frequency domain f and plot the single-sided amplitude spectrum P1. On average, longer signals produce better frequency approximations.
f = Fs*(0:(L/2))/L;

end