%% get_ABRthreshold.m
% 
% Estimate threshold amplitude of ABR response to tone burst that yields 
% 50% hit rate for detecting signal at input FPR / detection voltage RMS
%
% Searches space of Vd (detection ABR voltage RMS) using binary search
%
% Input: fpr - false positive rate (at zero ABR signal amplitude)
%        signal_V2   - ABR response RMS squared (m x 1 vector), derived from X where each of m examples is a column of length
%                           SAMPLES  (SAMPLES x m matrix) for ABR signal of
%                           zero amplitude
% Output: A_threshold - estimated ABR threshold (i.e. amplitude of ABR response 
%                       to signal input needed to yield 50% detection) 
%         Vd_est - estimated Vd detection voltage corresponding to input
%                  fpr
%
% Dependencies: simulateABR.m, get_ABRthreshold_Vd.m
% Last edit: 2/22/2019
%
% Author: George Liu

function [A_threshold, Vd_est] = get_ABRthreshold(fpr, signal_V2)

%%1) Get Vd from fpr
% search until estimated Vd gives FPR within this amt of input FPR 
PRECISION = 0.01;   
CLOSE_ENOUGH = 10^-8;

% initialize search parameters
Vd_min = 0;
Vd_max = 20;
error = 10^4;    % initialize FPR estimate to infinity

% find detection ABR voltage RMS (Vd) that corresponds to input FPR
while error > PRECISION
    Vd_mid = (Vd_min + Vd_max)/2;     % est ABR detection voltage^2 (uV^2)
    
    % Calculate hit rate
    hit = signal_V2>Vd_mid;
    miss = signal_V2<Vd_mid;
    hit_rate = sum(hit)/sum(hit+miss);
    
    % false positive rate is hit_rate for zero amplitude signal
    fpr_est = hit_rate;

	error = abs(fpr_est - fpr)/fpr;
    
    % update Vd_min or Vd_max
    if fpr_est - fpr > 0
        Vd_min = Vd_mid;
    else
        Vd_max = Vd_mid;
    end
    
    % end search if Vd_min and Vd_max have converged (because of duplicate 
    % signal_V2values
    if (Vd_max - Vd_min)/Vd_min < CLOSE_ENOUGH
        break
    end
end

Vd_est = Vd_mid;
% display(fpr_est)
% display(Vd_est)

%%2) Get A from Vd: Return ABR voltage amplitude threshold (uV) corresponding to Vd
% search until estimated Vd gives FPR within this amt of input FPR 
A_std = 1;
noise_mu = 0;
noise_sigma = 1;
m=1000;
A_threshold = get_ABRthreshold_Vd(Vd_est, A_std, noise_mu, noise_sigma, m);

end
