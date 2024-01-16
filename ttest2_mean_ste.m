function [tprob_2tail, tprob_1tail] = ttest2_mean_ste(n1, mean1, ste1, n2, mean2, ste2)
% Calculate 2 sample unpaired t-test for two samples (mean2 - mean1) with only their mean
% and standard errors.
%
% Input: 
%       n1 - sample size for group 1
%       mean1 - mean of group 1
%       ste1 - standard error of the estimated mean for group 1
%       n2 - sample size for group 2
%       mean2 - mean of group 2
%       ste2 - standard error of the estimated mean for group 2
%
% Output: 
%       tprob_2tail - 2 tailed probability for 2-sample unpaired t-test
%       tprob_1tail - 1 tailed probability for 2-sample unpaired t-test
%
% Dependencies: none
%
% Source: https://www.mathworks.com/matlabcentral/answers/301296-how-to-do-paired-t-test-with-mean-and-sd
%
% Last edit George Liu 9/15/2022

tval = (mean2 - mean1) / sqrt(ste1^2 + ste2^2);       % Calculate T-Statistic
v = (ste1^2 + ste2^2)^2 / (ste1^4/(n1 - 1) + ste2^4/(n2 - 1));      % Degrees of freedom
tdist2T = @(t, v) (1-betainc(v/(v+t^2),v/2,0.5));    % 2-tailed t-distribution
tdist1T = @(t, v) 1-(1-tdist2T(t,v))/2;              % 1-tailed t-distribution
tprob_2tail = 1-tdist2T(tval, v);
tprob_1tail = 1-tdist1T(tval, v);