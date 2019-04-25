% Generate fake ABR individual traces to check import_ABRcsv_folder.m
% 3-6-2019
% George Liu
% Code adapted from simulate_psychometric.m


% function [A, hit_rates] = simulate_psychometric(Vd, Amax, dA, A_std, noise_mu, noise_sigma, m)
Vd = 3;     % ABR detection voltage^2 (uV^2)
Amax = 9;  % signal amplitude (uV)
dA = 1;   % curve resolution
A_std = 1;
noise_mu = 0;
noise_sigma = 1;
m = 256;

A = (0:dA:Amax-dA)';    % vector of amplitudes (SPL) for simulating ABRs


%%Simulate ABR responses to different tone pip amplitudes to generate psychometric curve of
%%hit rates (TP+FP) vs SPL
n = numel(A);
X_many = cell(n,1);
hit_rates = zeros(n,1);
disp('Simulating psychometric function....')
for i = 1:n
    if mod(i,50*dA)==0
        disp(['Working on ', num2str(i/dA), ' of ', num2str(n/dA), '...'])
    end
%     hit_rates(i) = simulateABR(A(i), Vd, m);
    X_many{i} = simulateABR_many(A(i), A_std, noise_mu, noise_sigma, m); % (SAMPLES x m matrix)
    hit_rates(i) = sum(analyze_v2_ABR(X_many{i})>Vd)/m;
end