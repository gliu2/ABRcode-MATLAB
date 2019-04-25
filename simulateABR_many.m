%% simulateABR_many.m
% 
% Simulate distribution of ABR signals after repeated tone bursts
%
% Input: A_mu - average amplitude of ABR response
%        noise_A_mu - mean of normally distributed noise
%        noise_A_sigma - standard deviation of normally distributed noise
%
% Output: X - ABR responses, where each of m examples is a column of length
%               SAMPLES  (SAMPLES x m matrix)
%
% Dependencies: simulateABR_single.m
% Last edit: 2/7/2019
%
% Author: George Liu

function X = simulateABR_many(A_mu, A_std, noise_A_mu, noise_A_sigma, m)

% % A_mu = 4;   % average amplitude of ABR response
% A_std = 1;    % std ampltidue of ABR response
% noise_A_mu = 0;          % noise amplitude mean (uV)
% noise_A_sigma = 1;         % noise amplitude std dev (uV)
% m = 1000;    % number of ABR responses and tone pips played during experiment
SAMPLES = 1500;     % number of data points in one ABR response

% Collect distribution of <V^2>_t for signal ABR responses
X = zeros(SAMPLES, m);
for i = 1:m
    A = normrnd(A_mu, A_std);       % amplitude (uV)
    X(:,i) = simulateABR_single(A, noise_A_mu, noise_A_sigma);
end

end
