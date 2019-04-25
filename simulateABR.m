%% simulateABR.m
% 
% Simulate distribution of ABR signals after repeated tone bursts
%
% Input: A_mu - average amplitude of ABR response
%        Vd - detection threshold of ABR RMS squared (uV^2)
%
% Output: hit_rate - positive prediction probability given distribution of 
%                    simulated ABR RMS and input detection threshold (Vd)
%
% Dependencies: simulateABR_single.m, analyze_v2_ABR.m
% Last edit: 2/5/2019
%
% Author: George Liu

function hit_rate = simulateABR(A_mu, Vd, m)

% A_mu = 4;   % average amplitude of ABR response
A_std = 1;    % std ampltidue of ABR response
noise_A_mu = 0;          % noise amplitude mean (uV)
noise_A_sigma = 1;         % noise amplitude std dev (uV)

% Collect distribution of <V^2>_t for signal ABR responses
% m = 1000;    % iterations
signal_V2 = zeros(m, 1);
for i = 1:m
    A = normrnd(A_mu, A_std);       % amplitude (uV)
    x = simulateABR_single(A, noise_A_mu, noise_A_sigma);
    signal_V2(i,1) = analyze_v2_ABR(x);
end

% % Collect distribution of <V^2>_t for signal ABR responses
% n = 1000;    % iterations
% noise_V2 = zeros(n, 1);
% for i = 1:n
%     A = normrnd(0, 0);      % amplitude (uV)
%     x = simulateABR_single(A, noise_A_mu, noise_A_sigma);
%     noise_V2(i,1) = analyze_v2_ABR(x);
% end
% 
% combinedABR = [signal_V2; noise_V2];

% % plot histograms
% figure, histogram(signal_V2), xlabel('<V^2>_t (\muV)'), ylabel('Count'), title(['Signal ABR histogram (A=', num2str(A_mu), ')'])
% figure, histogram(noise_V2), xlabel('<V^2>_t (\muV)'), ylabel('Count'), title('Noise ABR histogram')
% figure, histogram(combinedABR), xlabel('<V^2>_t (\muV)'), ylabel('Count'), title(['Combined ABR histogram (A=', num2str(A_mu), ')'])

% xlabel('time (in seconds)');
% ylabel('ABR voltage (\muV)');
% title('Signal versus Time');

% Create psychometric-like function
% Vd = 3;   % ABR detection voltage^2
% TP = sum(signal_V2>Vd)
% FN = sum(signal_V2<Vd)
% FP = sum(noise_V2>Vd)
% TN = sum(noise_V2<Vd)
hit = signal_V2>Vd;
miss = signal_V2<Vd;
hit_rate = sum(hit)/sum(hit+miss);

end
