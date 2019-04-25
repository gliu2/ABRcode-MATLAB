%% plot_figure1.m
% 
% Plot of input sound signal
% Plot of output ABR response (signal, noise, combined, with histograms)
%
%
% Dependencies: simulateABR_many.m, get_ABRthreshold.m,
% simulate_psychometric.m
% Last edit: 2/25/2019
%
% Author: George Liu

rng(1) % seed random number generator to replicate results for debugging

% Get x axis of plots from simulateABR_single.m
%%Time specifications:
StopTime = 0.015;             % seconds
SAMPLES =1500;               % number of data points in single ABR response curve
Fs = SAMPLES/StopTime;          % samples per second
dt = 1/Fs;                   % seconds per sample
t = (0:dt:StopTime-dt)';     % seconds
latency = 0.003;            % seconds
Fc = 1000;                     % ABR frequencies: 500 Hz, 1 kHz, 2 kHz, 4 kHz

% input tone pip  
input_A_mu = 1;   % amplitude of ABR response
input_noise_A_mu = 0;          % noise amplitude mean (uV)
input_noise_A_sigma = 0;         % noise amplitude std dev (uV)
% X_in = simulateABR_single(input_A_mu, input_noise_A_mu, input_noise_A_sigma);
X_in = input_A_mu * sin(2*pi*Fc*(t-latency));

% Get experimental ABR data
A_mu = input_A_mu;   % average amplitude of ABR response
A_std = 1;    % std ampltidue of ABR response
noise_A_mu = 0;          % noise amplitude mean (uV)
noise_A_sigma = 1;         % noise amplitude std dev (uV)
m = 128;

X_signal = simulateABR_single(A_mu, 0, 0); % signal only
X_noise = simulateABR_single(0, noise_A_mu, noise_A_sigma); % noise only
X = X_signal + X_noise;  % Combined experimental ABR total output 

% For many repeated ABRs, get histograms
X_many = simulateABR_many(A_mu, A_std, noise_A_mu, noise_A_sigma, m); % SAMPLES x m matrix

TMIDDLE = 0.006; % s
TAPEX = 0.00474; % s
t_middle = find(t==TMIDDLE); % zero crossing near 0.006 ms for X_signal
t_apex = find(t==TAPEX); % negative peak near 0.0048 ms for X_signal

hist_middle = X_many(t_middle, :);
hist_apex = X_many(t_apex, :);

% Collect histograms of V RMS squared for many repeated ABRs
X_many_A0 = simulateABR_many(0, A_std, noise_A_mu, noise_A_sigma, m); % SAMPLES x m matrix

signal_V2 = analyze_v2_ABR(X_many);    % mx1 vector
noise_V2 = analyze_v2_ABR(X_many_A0);

% Plot psychometric curve
Vd = 3;     % ABR detection voltage^2 (uV^2)
Amax = 10;  % signal amplitude (uV)
dA = 0.2;   % curve resolution
[A, hit_rates] = simulate_psychometric(Vd, Amax, dA, A_std, noise_A_mu, noise_A_sigma, m);
FPR = hit_rates(1);

% Get Vd error for desired FPR
Nb = 10; % number bootstrap samples for est standard error of Vd, A_thresh
fpr_desired = FPR;
[A_thresh_std, Vd_std] = ABRthreshold_error(fpr_desired, X_many, Nb);

[A_plus1std, hit_rates_plus1std] = simulate_psychometric(Vd+Vd_std, Amax, dA, A_std, noise_A_mu, noise_A_sigma, m);
[A_minus1std, hit_rates_minus1std] = simulate_psychometric(Vd-Vd_std, Amax, dA, A_std, noise_A_mu, noise_A_sigma, m);
    

% Display key results
[A_thresh, Vd_est] = get_ABRthreshold(fpr_desired, X_many);
disp(['A_thresh: ', num2str(A_thresh), ' uV'])
disp(['A_thresh_std: ', num2str(A_thresh_std)])
disp(['Vd_est: ', num2str(Vd_est)])
disp(['Vd_std: ', num2str(Vd_std)])

%% Plot everything!
% Plot the input tone pip versus time:
figure;
plot(t,X_in);
xlabel('Time (seconds)', 'FontSize', 32);
ylabel('Amplitude (SPL)', 'FontSize', 32);
title(['Tone pip input (\omega=', num2str(Fc), ' Hz)'], 'FontSize', 40);
set(gca,'FontSize',20)
zoom xon;

% Plot the ABR response versus time:
figure;
plot(t,X_signal);
xlabel('Time (seconds)', 'FontSize', 32);
ylabel('ABR voltage (\muV)', 'FontSize', 32);
title('ABR response (signal only)', 'FontSize', 40);
set(gca,'FontSize',20)
zoom xon;

% Plot the ABR noise versus time:
figure;
plot(t,X_noise);
xlabel('Time (seconds)', 'FontSize', 32);
ylabel('ABR voltage (\muV)', 'FontSize', 32);
title('ABR noise', 'FontSize', 40);
set(gca,'FontSize',20)
zoom xon;

% Plot ABR response (signal and noise combined)
figure;
plot(t,X);
xlabel('Time (seconds)', 'FontSize', 32);
ylabel('ABR voltage (\muV)', 'FontSize', 32);
title('ABR response (signal and noise)', 'FontSize', 40);
set(gca,'FontSize',20)
zoom xon;

% % plot overlay of non-noisy ABR response
% hold on
% plot(t, x_nonoise, 'LineWidth', 1);
% hold off

% Histograms at middle and apex
figure, histogram(hist_middle)
xlabel('V (\muV)', 'FontSize', 32);
ylabel('Prob_V', 'FontSize', 32);
title(['ABR response at t=', num2str(TMIDDLE*1000), ' ms'], 'FontSize', 40)

figure, histogram(hist_apex)
xlabel('V (\muV)', 'FontSize', 32);
ylabel('Prob_V', 'FontSize', 32);
title(['ABR response at t=', num2str(TAPEX*1000), ' ms'], 'FontSize', 40)

% HIstograms of Voltage RMS squared
figure, histogram(signal_V2, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'FaceColor', 'r'), xlabel('<V^2>_t (\muV^2)', 'FontSize', 32), ylabel('Prob(<V^2>_t)', 'FontSize', 32), title('ABR histograms of RMS^2', 'FontSize', 40)
% figure, histogram(noise_V2, 'BinMethod', 'fd', 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 1, 'FaceAlpha', 0, 'FaceColor', 'b'), xlabel('<V^2>_t (\muV^2)', 'FontSize', 32), ylabel('Prob(<V^2>_t)', 'FontSize', 32), title('Noise ABR histogram', 'FontSize', 40)
hold on
histogram(noise_V2, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'FaceColor', 'b')
hold off
box off % remove ticks on top and right borders of plot
% axis tight % makes edges of data flush with left and right borders of plot
legend(['Signal (A=', num2str(A_mu), ')'],'Noise','location','northeast')
legend boxoff % remove box around legend

% Plot psychometric curve
%%Plot psychometric-like curve
figure, plot(A, hit_rates, 'LineWidth',2,'Color', 'blue')
% xlabel('SPL_{TH}')
xlabel('Signal ABR amplitude (\muV)', 'FontSize', 32)
ylabel('Probability of hit', 'FontSize', 32)
title(['Psychometric curve (Vd=', num2str(Vd), ' {\muV}^2, FP=', num2str(FPR), ')'], 'FontSize', 36)
hold on
plot(A_plus1std, hit_rates_plus1std, 'LineWidth',1,'Color', 'red')
plot(A_minus1std, hit_rates_minus1std, 'LineWidth',1,'Color', 'red')
hold off