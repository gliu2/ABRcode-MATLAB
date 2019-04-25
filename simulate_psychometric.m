%% simulate_psychometric.m
% 
% Create psychometric-like curve from simulated ABR responses to tone
% bursts
%
% Input: Vd - ABR detection voltage^2 (uV^2)
%        Amax - signal amplitude (uV)
%        dA - curve resolution
%        A_std -  standard deviation of ABR response amplitude
%        noise_A_mu - mean of normally distributed noise
%        noise_A_sigma - standard deviation of normally distributed noise
%        m - number of sweeps for each simulated ABR experiment; high number assures theoretical correct estimate
%
% Output: A - ABR signal amplitudes, 50x1 double array (x-axis)
%         hit_rates - detection probability given signal with ABR amplitude
%                     A, 50x1 double array (y-axis)
%
% Dependencies: simulateABR_many.m, analyze_v2_ABR.m
% Last edit: 2/22/2019
%
% Author: George Liu

function [A, hit_rates] = simulate_psychometric(Vd, Amax, dA, A_std, noise_mu, noise_sigma, m)
% Vd = 3;     % ABR detection voltage^2 (uV^2)
% Amax = 10;  % signal amplitude (uV)
% dA = 0.2;   % curve resolution

A = (0:dA:Amax-dA)';    % vector of amplitudes (SPL) for simulating ABRs


%%Simulate ABR responses to different tone pip amplitudes to generate psychometric curve of
%%hit rates (TP+FP) vs SPL
n = length(A);
hit_rates = zeros(n,1);
disp('Simulating psychometric function....')
for i = 1:n
    if mod(i,50*dA)==0
        disp(['Working on ', num2str(i/dA), ' of ', num2str(n/dA), '...'])
    end
%     hit_rates(i) = simulateABR(A(i), Vd, m);
    X = simulateABR_many(A(i), A_std, noise_mu, noise_sigma, m);
    hit_rates(i) = sum(analyze_v2_ABR(X)>Vd)/m;
end

% % Calculate FPR unique to Vd
% FPR = hit_rates(1);
% 
% %%Plot psychometric-like curve
% figure, plot(A, hit_rates)
% % xlabel('SPL_{TH}')
% xlabel('Signal ABR amplitude (\muV)')
% ylabel('Probability of hit')
% title(['Psychometric curve for signal ABR detection (Vd=', num2str(Vd), ' {\muV}^2, FP=', num2str(FPR), ')'])

end
