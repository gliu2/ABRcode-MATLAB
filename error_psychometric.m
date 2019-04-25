%% error_psychometric.m
% 
% Create psychometric-like curve from simulated ABR responses to tone
% bursts
%
% Input: Vd - ABR detection voltage^2 (uV^2)
%        A - desired ABR response amplitude for estimating psychometric
%            P(hit) error
%        A_std -  standard deviation of ABR response amplitude
%        noise_A_mu - mean of normally distributed noise
%        noise_A_sigma - standard deviation of normally distributed noise
%        m - number of sweeps for each simulated ABR experiment; high number assures theoretical correct estimate
%        nsample - number of samples to estimate P(hit) std
%
% Output: hit_rates_std - std of detection probability given signal with ABR amplitude
%                         A, double array of size length(A)x1
%
% Dependencies: simulateABR_many.m, analyze_v2_ABR.m
% Last edit: 2/27/2019
%
% Author: George Liu

function hit_rates_std = error_psychometric(Vd, A, A_std, noise_mu, noise_sigma, m, nsample)
% Vd = 3;     % ABR detection voltage^2 (uV^2)
% Amax = 10;  % signal amplitude (uV)
% dA = 0.2;   % curve resolution

%%Simulate ABR responses to different tone pip amplitudes to generate psychometric curve of
%%hit rates (TP+FP) vs SPL
n = length(A);
hit_rates_std = zeros(n, 1);
for i = 1:n
    hit_rates = zeros(nsample,1);
    disp(['Estimating psychometric function error at A=', num2str(A(i)), ' ....'])
    for j = 1:nsample
        if mod(j,100)==0
            disp(['Working on ', num2str(j), ' of ', num2str(nsample), '...'])
        end
    %     hit_rates(i) = simulateABR(A(i), Vd, m);
        X = simulateABR_many(A(i), A_std, noise_mu, noise_sigma, m);
        hit_rates(j) = sum(analyze_v2_ABR(X)>Vd)/m;
    end
    hit_rates_std(i) = std(hit_rates);
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