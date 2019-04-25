%% get_ABRthreshold_Vd.m
% 
% Estimate theoretical ABR threshold of signal intensity necessary to
% produce ABR "hit" 50% of time.
%
% Input: Vd - RMS^2 detection threshold (sets FPR at zero ABR signal amplitude)
%        A_std -  standard deviation of ABR response amplitude
%        noise_A_mu - mean of normally distributed noise
%        noise_A_sigma - standard deviation of normally distributed noise
%        m - number of sweeps for each simulated ABR experiment; high number assures theoretical correct estimate
%
% Output: A_threshold - estimated ABR threshold (i.e. amplitude of ABR response 
%                       to signal input needed to yield 50% detection) 
%         A_cache - cached amplitudes to plot psychometric curve (?x1 vector)
%         hitrate_cache - cached hit-rates to plot psychometric curve (?x1 vector)
%
% Dependencies: simulateABR_many.m, analyze_v2_ABR.m
% Last edit: 2/22/2019
%
% Author: George Liu

function [A_threshold, A_cache, hitrate_cache] = get_ABRthreshold_Vd(Vd, A_std, noise_mu, noise_sigma, m)

%%Get A from Vd: Return ABR voltage amplitude threshold (uV) corresponding to Vd
% search until estimated Vd gives FPR within this amt of input FPR 
PRECISION2 = 0.001;   
TARGET_HITRATE = 0.5;
CLOSE_ENOUGH = 10^-8; % end after many futile iterations that bring binary search box less than this width; deals with overlapping values for bootstrapping
% m=1000; % number of sweeps for each ABR experiment; high number assures theoretical correct estimate

% initialize search parameters
A_min = 0;
A_max = 20;
error2 = 10^4;    % initialize A_mid estimate to infinity

% find ABR amplitude that yields 50% pred at Vd using binary search
A_cache = [];
hitrate_cache = [];
while error2 > PRECISION2
    A_mid = (A_max + A_min)/2;
    % estimate hitrate for A_mid, Vd
    X = simulateABR_many(A_mid, A_std, noise_mu, noise_sigma, m);
    est_hitrate = sum(analyze_v2_ABR(X)>Vd)/m;
    
    % compute error
	error2 = abs(est_hitrate - TARGET_HITRATE)/TARGET_HITRATE;
    
%     disp(['Amax: ', num2str(A_max), ', Amin: ', num2str(A_min), ', Amid: ', num2str(A_mid)]) % debugging
    
    % update A_min or A_max
    if est_hitrate > TARGET_HITRATE
        A_max = A_mid;
    else
        A_min = A_mid;
    end
    
    % Cache searched amplitudes and hit rates to plot psychometric curve
    A_cache = [A_cache; A_mid];
    hitrate_cache = [hitrate_cache; est_hitrate];
    
    % end search if Vd_min and Vd_max have converged (because of duplicate 
    % signal_V2values
    if (A_max - A_min)/(A_min + CLOSE_ENOUGH) < CLOSE_ENOUGH
        break
    end
end

% Output best A
A_threshold = A_mid;

end