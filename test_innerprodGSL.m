% test inner product GSL
% George Liu
% 5-30-19
% Dependencies:

SAMPLES = size(X_csv{1}, 1);
SAMPLING_RATE = 200000; % Hz
dt = 1/SAMPLING_RATE * 1000; % ms
T_total = (SAMPLES-1)*dt; % ms
x = 0:dt:T_total;

% %% Obtain inner product distributions with Fourier transforms
% [f, fourier_innerproddist] = fourier_innerprod_ABR(X_csv, SAMPLING_RATE);
% [p_val, ranksum_stat, p_val_KS, ks_stat] = innerprod2pval(fourier_innerproddist);
% visualize_innerprod_stats(p_val, ranksum_stat, p_val_KS, ks_stat, A_csv)
% visualize_innerprod_hist(X_csv, x, A_csv, fourier_innerproddist)

%% Obtain inner product distributions in time domain
% chunk = 
avg_lag_xcovExtrap = lags_xcov(X_csv, dt, 0.5, A_csv, true, chunk);