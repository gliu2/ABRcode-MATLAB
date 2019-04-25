%% get_ABRthreshold_se.m
% 
% Estimate standard error of ABR detection threshold (Vd) estimate using
% synthetic ABR experiment data for an input FPR
%
% Uses boot strapping
%
% Input: fpr - false positive rate (at zero ABR signal amplitude)
% Output: A_threshold - estimated ABR threshold (i.e. amplitude of ABR response 
%                       to signal input needed to yield 50% detection) 
%
% Dependencies: simulateABR_many.m, get_ABRthreshold.m
% Last edit: 2/7/2019
%
% Author: George Liu

% function A_threshold = get_ABRthreshold_se(fpr)

rng(1) % seed random number generator to replicate results for debugging

% Get experimental ABR data
A_mu = 0;   % average amplitude of ABR response
A_std = 1;    % std ampltidue of ABR response
noise_A_mu = 0;          % noise amplitude mean (uV)
noise_A_sigma = 1;         % noise amplitude std dev (uV)
m = 1000;

X = simulateABR_many(A_mu, A_std, noise_A_mu, noise_A_sigma, m); % SAMPLES x m matrix
signal_V2 = analyze_v2_ABR(X);    % mx1 vector

% bootstrap to draw Nb samples with replacement
Nb = 100;
boot_indexes = ceil(rand(m, Nb)*m);
boot_samples = signal_V2(boot_indexes); % m x Nb matrix, columns are bootstrap samples

% Estimate Vd for each boot strap sample
fpr = 0.003;
Vd_est = zeros(Nb, 1);
A_threshold = zeros(Nb, 1);
for i=1:Nb
    disp(['Working on ', num2str(i), ' out of ', num2str(Nb), ' bootstrap samples....'])
    this_bootsample = boot_samples(:,i);
    [A_threshold(i,1), Vd_est(i,1)] = get_ABRthreshold(fpr, this_bootsample);
end    

% Calculated standard error of Vd statistic from boot strapping
Vd_mean = mean(Vd_est)
Vd_std = sqrt(sum((Vd_est - Vd_mean).^2)/(Nb-1))

% Experimental Vd statistic 
[A_threshold_exp, Vd_est_exp] = get_ABRthreshold(fpr, signal_V2);
display(Vd_est_exp)

% end
