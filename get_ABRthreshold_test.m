%% get_ABRthreshold_test.m
% 
% Estimate threshold amplitude of ABR response to tone burst that yields 
% 50% hit rate for detecting signal at input FPR / detection voltage RMS
%
% Searches space of Vd (detection ABR voltage RMS) using binary search
%
% Input: fpr - false positive rate (at zero ABR signal amplitude)
% Output: A_threshold - estimated ABR threshold (i.e. amplitude of ABR response 
%                       to signal input needed to yield 50% detection) 
%
% Dependencies: simulateABR_many.m, analyze_v2_ABR.m
% Last edit: 2/7/2019
%
% Author: George Liu

function A_threshold = get_ABRthreshold_test(fpr)

% Get experimental ABR data
A_mu = 0;   % average amplitude of ABR response
A_std = 1;    % std ampltidue of ABR response
noise_A_mu = 0;          % noise amplitude mean (uV)
noise_A_sigma = 1;         % noise amplitude std dev (uV)
m = 1000;

X = simulateABR_many(A_mu, A_std, noise_A_mu, noise_A_sigma, m); % SAMPLES x m matrix
signal_V2 = analyze_v2_ABR(X);    % mx1 vector

[A_threshold, Vd_est] = get_ABRthreshold(fpr, signal_V2);

display(Vd_est)
display(A_threshold)


% %%1) Get Vd from fpr
% % search until estimated Vd gives FPR within this amt of input FPR 
% PRECISION = 0.01;   
% 
% % initialize search parameters
% Vd_min = 0;
% Vd_max = 20;
% error = 10^4;    % initialize FPR estimate to infinity
% 
% % find detection ABR voltage RMS (Vd) that corresponds to input FPR
% while error > PRECISION
%     Vd_mid = (Vd_min + Vd_max)/2;     % est ABR detection voltage^2 (uV^2)
%     fpr_est = simulateABR(0, Vd_mid);
% 	error = abs(fpr_est - fpr)/fpr;
%     
%     % update Vd_min or Vd_max
%     if fpr_est - fpr > 0
%         Vd_min = Vd_mid;
%     else
%         Vd_max = Vd_mid;
%     end
% end
% 
% Vd_est = Vd_mid;
% display(fpr_est)
% display(Vd_est)
% 
% %%2) Get A from Vd: Return ABR voltage amplitude threshold (uV) corresponding to Vd
% % search until estimated Vd gives FPR within this amt of input FPR 
% PRECISION2 = 0.01;   
% TARGET_HITRATE = 0.5;
% 
% % initialize search parameters
% A_min = 0;
% A_max = 20;
% error2 = 10^4;    % initialize A_mid estimate to infinity
% 
% % find ABR amplitude that yields 50% pred at Vd using binary search
% while error2 > PRECISION2
%     A_mid = (A_max + A_min)/2;
%     % estimate hitrate for A_mid, Vd
%     est_hitrate = simulateABR(A_mid, Vd_est);
%     
%     % compute error
% 	error2 = abs(est_hitrate - TARGET_HITRATE)/TARGET_HITRATE;
%     
%     % update A_min or A_max
%     if est_hitrate > TARGET_HITRATE
%         A_max = A_mid;
%     else
%         A_min = A_mid;
%     end
% end
% 
% % Output best A
% A_threshold = A_mid;

end
