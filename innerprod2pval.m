%% innerprod2pval.m
% 
% Wilcoxin sign rank and Kolmogorov-Smirnov two-sample test for input
% inner product distributions. Tests if nonzero dB distributions are 
% different from 0 dB single trace distribution. Returns p-value per 
% stimulus level to estimate ABR threshold.
%
% Input: dist_innerprod - Inner product distributions (lag-adjusted).
%               Cell of size (A_length, 1), where A_length is # of dB levels.
%               Entry dist_innerprod{i} gives vector of size (m_traces, 1),
%               whose entry j is inner product (lag-adjusted) of j-th
%               single trace (at i-th dB level) with max dB level average ABR trace.
%
% Output: p_val - p-value from Wilcoxin sign rank two-sample test, 2-sided
%                 Vector of size (A_length, 1).
%
%         signrank_stat - test statistic from Wilcoxin sign rank test
%                         Vector of size (A_length, 1).
%
%         p_val_KS - p-value from Kolmogorov?Smirnov two-sample test
%                    Vector of size (A_length, 1).
%
%         ks_stat - test statistic from Kolmogorov?Smirnov test
%                   Vector of size (A_length, 1).
%
% Dependencies: none
% Last edit: 5/26/2019
%
% Author: George Liu

function [p_val, signrank_stat, p_val_KS, ks_stat] = innerprod2pval(dist_innerprod)
A_length = size(dist_innerprod, 1);

zero_dB_innerdist = dist_innerprod{1};
p_val = zeros(A_length, 1);
signrank_stat = zeros(A_length, 1);
p_val_KS = zeros(A_length, 1);
ks_stat = zeros(A_length,1);
for i = 1:A_length
    [p_val(i), ~, stats] = signrank(dist_innerprod{i}, zero_dB_innerdist); % Wilcoxin sign-rank test, two-sided 
    signrank_stat(i) = stats.signedrank;
    [~, p_val_KS(i), ks_stat(i)] = kstest2(dist_innerprod{i}, zero_dB_innerdist); % Two-sample Kolmogorov-Smirnov test (unequal cdf's)
end

end