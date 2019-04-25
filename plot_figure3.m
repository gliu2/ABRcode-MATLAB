%% plot_figure1.m
% 
% Plot histograms of large sample limit w/ and w/o averaging
%
% Dependencies: simulateABR_many.m, get_ABRthreshold.m,
% simulate_psychometric.m
% Last edit: 2/25/2019
%
% Author: George Liu

rng(1) % seed random number generator to replicate results for debugging

% ABR response parameters
A_mu = 1;   % average amplitude of ABR response
A_std = 1;    % std ampltidue of ABR response
A_std_noise = 0; % std ampltidue of ABR response for calculating noise
noise_A_mu = 0;          % noise amplitude mean (uV)
noise_A_sigma = 10;         % noise amplitude std dev (uV)
m = 128; % experimental
m_limit = 3000; % large sample limit

% Simulate repeated ABR responses - single traces
X_many_exp = simulateABR_many(A_mu, A_std, noise_A_mu, noise_A_sigma, m); % SAMPLES x m matrix
X_many_largem = simulateABR_many(A_mu, A_std, noise_A_mu, noise_A_sigma, m_limit); % SAMPLES x m_limit matrix

% X_many_A0 = simulateABR_many(0, 0, noise_A_mu, noise_A_sigma, m); % SAMPLES x m matrix
X_many_A0 = simulateABR_many(0, A_std_noise, noise_A_mu, noise_A_sigma, m); % SAMPLES x m matrix
% X_many_A0_largem = simulateABR_many(0, 0, noise_A_mu, noise_A_sigma, m_limit); % SAMPLES x m_limit matrix
X_many_A0_largem = simulateABR_many(0, A_std_noise, noise_A_mu, noise_A_sigma, m_limit); % SAMPLES x m_limit matrix

% Simulate repeated ABR responses - averaged trace
Vrms2_exp_ave = zeros(m, 1);
Vrms2_A0_ave = zeros(m, 1);
for i=1:m
    X_exp_i = simulateABR_many(A_mu, A_std, noise_A_mu, noise_A_sigma, m); % SAMPLES x m matrix
    X_exp_ave = mean(X_exp_i, 2); % SAMPLES x 1 matrix
    Vrms2_exp_ave(i) = analyze_v2_ABR(X_exp_ave); % 1x1 vector
    
    X_A0_i = simulateABR_many(0, A_std_noise, noise_A_mu, noise_A_sigma, m); % SAMPLES x m matrix
    X_A0_ave = mean(X_A0_i, 2); % SAMPLES x 1 matrix
    Vrms2_A0_ave(i) = analyze_v2_ABR(X_A0_ave); % 1x1 vector
    
end

% X_many_exp_ave = mean(X_many_exp, 2); % SAMPLES x 1 matrix
% X_many_largem_ave = mean(X_many_largem, 2); % SAMPLES x 1 matrix
% X_many_A0_ave = mean(X_many_A0, 2); % SAMPLES x 1 matrix
% X_many_A0__largem_ave = mean(X_many_A0_largem, 2); % SAMPLES x 1 matrix

% RMS squared
Vrms2_exp_single = analyze_v2_ABR(X_many_exp);    % mx1 vector
Vrms2_largem_single = analyze_v2_ABR(X_many_largem); % m_limit x 1 vector
Vrms2_A0_single = analyze_v2_ABR(X_many_A0); % mx1 vector
Vrms2_A0_largem_single = analyze_v2_ABR(X_many_A0_largem); % m_limit x 1 vector


%% Plot
% HIstograms of Voltage RMS squared
figure, histogram(Vrms2_exp_single, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'FaceColor', 'r')
xlabel('{V_{RMS}}^2 (\muV^2)', 'FontSize', 32)
ylabel('Prob({V_{RMS}}^2)', 'FontSize', 32)
title(['ABR histograms of {V_{RMS}}^2 (', num2str(m), ' single traces)'], 'FontSize', 40)
% figure, histogram(noise_V2, 'BinMethod', 'fd', 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 1, 'FaceAlpha', 0, 'FaceColor', 'b'), xlabel('<V^2>_t (\muV^2)', 'FontSize', 32), ylabel('Prob(<V^2>_t)', 'FontSize', 32), title('Noise ABR histogram', 'FontSize', 40)
hold on
histogram(Vrms2_A0_single, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'FaceColor', 'b')
hold off
box off % remove ticks on top and right borders of plot
% axis tight % makes edges of data flush with left and right borders of plot
legend(['Signal (A=', num2str(A_mu), '), m=', num2str(1)], ['Noise, m=', num2str(1)], 'location', 'northeast')
legend boxoff % remove box around legend

% HIstograms of Voltage RMS squared
figure, histogram(Vrms2_largem_single, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'FaceColor', 'r')
xlabel('{V_{RMS}}^2 (\muV^2)', 'FontSize', 32)
ylabel('Prob({V_{RMS}}^2)', 'FontSize', 32)
title(['ABR histograms of {V_{RMS}}^2 (', num2str(m_limit), ' single traces)'], 'FontSize', 40)
% figure, histogram(noise_V2, 'BinMethod', 'fd', 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 1, 'FaceAlpha', 0, 'FaceColor', 'b'), xlabel('<V^2>_t (\muV^2)', 'FontSize', 32), ylabel('Prob(<V^2>_t)', 'FontSize', 32), title('Noise ABR histogram', 'FontSize', 40)
hold on
histogram(Vrms2_A0_largem_single, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'FaceColor', 'b')
hold off
box off % remove ticks on top and right borders of plot
% axis tight % makes edges of data flush with left and right borders of plot
legend(['Signal (A=', num2str(A_mu), '), m=', num2str(1)], ['Noise, m=', num2str(1)], 'location', 'northeast')
legend boxoff % remove box around legend

% HIstograms of Voltage RMS squared
figure, histogram(Vrms2_exp_ave, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'FaceColor', 'r')
xlabel('{V_{RMS}}^2 (\muV^2)', 'FontSize', 32)
ylabel('Prob({V_{RMS}}^2)', 'FontSize', 32)
title('ABR histograms of {V_{RMS}}^2, old method (averaged traces)', 'FontSize', 40)
% figure, histogram(noise_V2, 'BinMethod', 'fd', 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 1, 'FaceAlpha', 0, 'FaceColor', 'b'), xlabel('<V^2>_t (\muV^2)', 'FontSize', 32), ylabel('Prob(<V^2>_t)', 'FontSize', 32), title('Noise ABR histogram', 'FontSize', 40)
hold on
histogram(Vrms2_A0_ave, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'FaceColor', 'b')
hold off
box off % remove ticks on top and right borders of plot
% axis tight % makes edges of data flush with left and right borders of plot
legend(['Signal (A=', num2str(A_mu), '), m=', num2str(m)], ['Noise, m=', num2str(m)], 'location', 'northeast')
legend boxoff % remove box around legend