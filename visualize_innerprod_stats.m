%% visualize_innerprod_stats.m
% 
% Plot Wilcoxin sign rank and Kolmogorov-Smirnov two-sample test p-values
% and test statistics. Can be run on output of innerprod2pval.m.
%
% Input: p_val - p-value from Wilcoxin sign rank two-sample test, 2-sided
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
%         A_csv - stimulus levels (db SPL).
%                 Vector of size (A_length, 1).
%
% Output: 1 Figure with 4 subplots of p-values and test statistics per dB.
%
% Dependencies: none
% Last edit: 5/26/2019
%
% Author: George Liu

function visualize_innerprod_stats(p_val, signrank_stat, p_val_KS, ks_stats, A_csv)

THRESH_INDEX = 5; % index of A_csv value that yields p_val < 0.05; from looking at plot

figure('DefaultAxesFontSize', 20)
% Plot p-value vs dB SPL level: Wilcoxin sign rank test
subplot(2,2,1)
plot(A_csv, p_val, '-o')
title(['Wilcoxin sign rank for inner products vs 0 dB SPL'])
xlabel('Amplitude (dB SPL)')
ylabel('p-value')  
% CUSTOM add coordinate of threshold point
x_text = A_csv(THRESH_INDEX);
y_text = p_val(THRESH_INDEX);
text(x_text - 8, y_text, ['(' num2str(x_text) ', ' num2str(y_text, 2) ')'])

% Plot Wilcoxign sign rank z-stat vs dB SPL level: Wilcoxin sign rank test
subplot(2,2,2)
plot(A_csv, signrank_stat, '-o')
title(['Wilcoxin sign rank for inner products vs 0 dB SPL'])
xlabel('Amplitude (dB SPL)')
ylabel('Sign rank test statistic')  
% CUSTOM add coordinate of threshold point
THRESH_INDEX = 5; % index of A_csv value that yields p_val < 0.05; from looking at plot
x_text = A_csv(THRESH_INDEX);
y_text = signrank_stat(THRESH_INDEX);
text(x_text - 8, y_text, ['(' num2str(x_text) ', ' num2str(y_text, 2) ')'])

% Plot p-value vs dB SPL level: K-S test
subplot(2,2,3)
plot(A_csv, p_val_KS, '-o')
title(['Kolmogorov-Smirnov test for inner products vs 0 dB SPL'])
xlabel('Amplitude (dB SPL)')
ylabel('p-value')  
% CUSTOM add coordinate of threshold point
THRESH_INDEX = 5; % index of A_csv value that yields p_val < 0.05; from looking at plot
x_text = A_csv(THRESH_INDEX);
y_text = p_val_KS(THRESH_INDEX);
text(x_text - 8, y_text, ['(' num2str(x_text) ', ' num2str(y_text, 2) ')'])

% Plot KS_score vs dB SPL level: K-S test
subplot(2,2,4)
plot(A_csv, ks_stats, '-o')
title(['Kolmogorov-Smirnov test for inner products vs 0 dB SPL'])
xlabel('Amplitude (dB SPL)')
ylabel('K-S score')  
% CUSTOM add coordinate of threshold point
THRESH_INDEX = 5; % index of A_csv value that yields p_val < 0.05; from looking at plot
x_text = A_csv(THRESH_INDEX);
y_text = ks_stats(THRESH_INDEX);
text(x_text - 8, y_text, ['(' num2str(x_text) ', ' num2str(y_text, 2) ')'])

end