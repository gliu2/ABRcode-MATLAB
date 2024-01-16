% Script to implement blackman filter to signal
FONTSIZE = 14;
M = 64; % window length
N = 250;

% initialize fake signal
% Signal parameters (synthetic sinusoid):
wxT = 2*pi/N;% Sinusoid frequency in rad/sample
A = 1;          % Sinusoid amplitude
phix = pi/3;    % Sinusoid phase

% Compute the signal x:
n = 0:N-1;    % time indices for sinusoid and FFT
x = A * cos(wxT*n+phix); % a sinusoid

% Implement filter
win = blackman(M);
win = win/sum(win); % Normalization
% wvtool(win)
filtered = conv(x,win,'same');

fig = figure;
t = tiledlayout(1, 2, 'TileIndexing', 'rowmajor');
nexttile(t)
plot(n, x, 'blue', 'LineWidth', 3)
xlabel('Index')
ylabel('Amplitude')
title('Synthetic signal')
set(gca,'FontSize', FONTSIZE)

nexttile(t)
plot(1:numel(filtered), filtered, 'blue', 'LineWidth', 3)
xlabel('Index')
ylabel('Amplitude')
title('Blackman filtered signal')
set(gca,'FontSize', FONTSIZE)