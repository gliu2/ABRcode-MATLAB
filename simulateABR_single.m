%% simulateABR_single.m
% 
% Simulate ABR signal for a single tone burst
%
% ABR parameters obtained from: https://www.audiologyonline.com/articles/update-on-auditory-evoked-responses-19040
%
% Input: A - amplitude of ABR response sine wave 
%        noise_A_mu - mean of normally distributed noise
%        noise_A_sigma - standard deviation of normally distributed noise
%
% Output: x - ABR voltage output per time increment dt (nx1 vector)
%
% Dependencies: none
% Last edit: 2/6/2019
%
% Author: George Liu

function x = simulateABR_single(A, noise_A_mu, noise_A_sigma)

%%Time specifications:
StopTime = 0.015;             % seconds
SAMPLES =1500;               % number of data points in single ABR response curve
Fs = SAMPLES/StopTime;          % samples per second
dt = 1/Fs;                   % seconds per sample
t = (0:dt:StopTime-dt)';     % seconds
latency = 0.003;            % seconds

%%Other wave parameters
% A = 1;                       % amplitude (uV)
% noise_A_mu = 0;          % noise amplitude mean (uV)
% noise_A_sigma = 1;         % noise amplitude std dev (uV)

%%Sine wave:
Fc = 1000;                     % ABR frequencies: 500 Hz, 1 kHz, 2 kHz, 4 kHz
x_nolatency = A*sin(2*pi*Fc*(t-latency));

%%simulate ABR response to tone burst
filter = heaviside(t-latency);
x_nonoise = x_nolatency.*filter;
noise = normrnd(noise_A_mu, noise_A_sigma, size(t));
x = x_nonoise + noise;


% % Plot the signal versus time:
% figure;
% plot(t,x);
% xlabel('time (in seconds)');
% ylabel('ABR voltage (\muV)');
% title('Signal versus Time');
% zoom xon;
% 
% % plot overlay of non-noisy ABR response
% hold on
% plot(t, x_nonoise, 'LineWidth', 1);
% hold off

end
