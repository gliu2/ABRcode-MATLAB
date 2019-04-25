%% plot_figure2.m
% 
% Plot of Thresh_new vs Thresh_old
% Compare sensitivity of ABR measurement with new and old thresholds
%
% Dependencies: simulateABR_many.m, get_ABRthreshold.m,
% simulate_psychometric.m
% Last edit: 2/11/2019
%
% Author: George Liu

rng(2) % seed random number generator to replicate results for debugging

% Get experimental ABR data
A_mu = 1;   % average amplitude of ABR response
A_std = 1;    % std ampltidue of ABR response
noise_mu = 0;          % noise amplitude mean (uV)
noise_sigma = 1;         % noise amplitude std dev (uV)
% m = 128;
% mrange = 128:128:1000;
mrange = 128;
err_Ath_old_cache = zeros(numel(mrange), 1);
err_Ath_new_cache = zeros(numel(mrange), 1);
Vd_cache = zeros(numel(mrange), 1);
Vd_old_cache = zeros(numel(mrange), 1);
V_0_cache = zeros(numel(mrange), 1);
count = 0;
for m = mrange
    count = count + 1;
    disp(['m: ', num2str(m)])

    % For many repeated ABRs, get histograms
    X_many = simulateABR_many(A_mu, A_std, noise_mu, noise_sigma, m); % SAMPLES x m matrix
    X_many_noise = simulateABR_many(0,0, noise_mu, noise_sigma, m); % SAMPLES x m matrix

    % TMIDDLE = 0.006; % s
    % TAPEX = 0.00474; % s
    % t_middle = find(t==TMIDDLE); % zero crossing near 0.006 ms for X_signal
    % t_apex = find(t==TAPEX); % negative peak near 0.0048 ms for X_signal
    % 
    % hist_middle = X_many(t_middle, :);
    % hist_apex = X_many(t_apex, :);    
    % 
    % % Collect histograms of V RMS squared for many repeated ABRs
    % X_many_A0 = simulateABR_many(0, A_std, noise_A_mu, noise_A_sigma, m); % SAMPLES x m matrix


    % Old threshold shoots for FPR 1% or Fsp 3.1, SNR 

    % Simulate averaged V for synthetic Ricci lab experimental data
    X_avg = mean(X_many, 2); % SAMPLES x 1 vector of averaged ABR responses
    X_noise_avg = mean(X_many_noise, 2); % SAMPLES x 1 vector of averaged ABR response to noise

    V_0 = sqrt(analyze_v2_ABR(X_noise_avg)); % RMS no signal; also standard deviation of noise averaged
    V_SPL = sqrt(analyze_v2_ABR(X_avg)); % RMS signal

    % sigma_noise_avg = std(X_noise_avg); % RMS no signal std

    SNR = V_SPL / V_0; % assumimng uncorrelated noise
    Fsp = std(X_avg).^2/V_0.^2; % uses individual trace data

    MULT = 2; % # std of noise
    Vd_RMS = (MULT*V_0)*sqrt(m); % RMS threshold based on old estimate (sets FPR for old and new method)
    Vd = Vd_RMS.^2; % new detection threshold for RMS squared
    Vd_old = (MULT*V_0).^2; % old detection threshold for RMS squared
    disp(Vd)
    disp(SNR)
    disp(Fsp)

    % Get old ABR threshold
    Amax = 10;
    dA = 0.2;
    A = (0:dA:Amax-dA)';    % vector of amplitudes (SPL) for simulating ABRs
    n = length(A);
    hit_rates_old = zeros(n,1);
    hit_rates_new = zeros(n,1);
    V_i = zeros(m,n);
    disp('Simulating psychometric function....')
    for i = 1:n
        if mod(i,50*dA)==0
            disp(['Working on ', num2str(i/dA), ' of ', num2str(n/dA), '...'])
        end
        X_i = simulateABR_many(A(i),A_std, noise_mu, noise_sigma, m); % SAMPLES x m matrix
        X_i_avg = mean(X_i, 2); % SAMPLES x 1 vector of averaged ABR response to noise
        V_i_avg = analyze_v2_ABR(X_i_avg);
        V_i(:,i) = analyze_v2_ABR(X_i);
    %     disp(V_i(1:5, i))
        hit_rates_old(i) = sum(V_i_avg > Vd_old);
        hit_rates_new(i) = sum(V_i(:, i) > Vd)/m;
    %     disp(analyze_v2_ABR(X_i))
    %     disp(hit_rates_new(i))
    end

    % Get theoretical ABR psychometric curve and threshold
    m_theor = 3000; % number of samples to generate theoretical psychometric curve
    [A_theor, hitrates_theor] = simulate_psychometric(Vd, Amax, 0.2, A_std, noise_mu, noise_sigma, m_theor);
    disp('Calculate get ABR threshold....')
    [A_threshold_theor, A_cache, hitrate_cache] = get_ABRthreshold_Vd(Vd, A_std, noise_mu, noise_sigma, m_theor);
    disp('Done calculating.')
    A_theoretical = [A_theor; A_cache];
    hitrate_theoretical = [hitrates_theor; hitrate_cache];

    % put x and y coords in order of ascending x
    [A_theoretical,I] = sort(A_theoretical);
    hitrate_theoretical = hitrate_theoretical(I);


    % Get "experimental" ABR threshold
    A_exp_old = ZeroGL(A, hit_rates_old-0.5);  % not sure why also returns  9.7000
    if isempty(A_exp_old)
        A_exp_old = A(end); % if no hits, set threshold to maximum possible value
    else 
        A_exp_old = A_exp_old(1);
    end
    A_exp_new = ZeroGL(A, hit_rates_new-0.5);  % not sure why also returns  9.7000
    if isempty(A_exp_new)
        A_exp_new = A(end); % if no hits, set threshold to maximum possible value
    else 
        A_exp_new = A_exp_new(1);
    end
    all_thresholds = [A_exp_old, A_exp_new, A_threshold_theor]; % for plotting 

    disp(['A threshold theoretical: ', num2str(A_threshold_theor)])
    
    % Calculate error bars for psychometric curves
    nsample = 1000;
    A_errorbars = 0.5:0.5:9.5;
    disp('Calculate error bars for psychometric curves.')
    hit_rates_std_exp = error_psychometric(Vd, A_errorbars, A_std, noise_mu, noise_sigma, m, nsample);
    hit_rates_std_theor = error_psychometric(Vd, A_errorbars, A_std, noise_mu, noise_sigma, m_theor, nsample);
    hit_rates_std_old = error_psychometric(Vd_old, A_errorbars, A_std, noise_mu, noise_sigma, m, nsample);
    
    k = round(A_errorbars/dA+1, 0);
    k_theor = zeros(length(A_errorbars), 1);
    error_bar_exp = NaN(size(A));
    error_bar_theor = NaN(size(A_theoretical));
    error_bar_old = NaN(size(A));
    for j=1:length(A_errorbars)
        disp(num2str(j))
        k_theor(j) = round(find(round(A_theoretical, 1)==round(A_errorbars(j), 1), 1), 0);
        disp('c')
        error_bar_exp(k(j)) = hit_rates_std_exp(j);
        error_bar_theor(k_theor(j)) = hit_rates_std_theor(j);
        error_bar_old(k(j)) = hit_rates_std_old(j);
    end

    % Calculate error of estimate
    err_Ath_old = abs(A_threshold_theor - A_exp_old)/A_threshold_theor;
    err_Ath_new = abs(A_threshold_theor - A_exp_new)/A_threshold_theor;
    disp(['Ath old rel error: ', num2str(err_Ath_old)])
    disp(['Ath new rel error: ', num2str(err_Ath_new)])

%     %%Plot psychometric-like curve
%     h=figure; plot(A, hit_rates_new, 'LineWidth', 1, 'Color', 'red')
%     hold on
%     plot(A, hit_rates_old, 'LineWidth', 1, 'Color', 'blue')
%     plot(A_theoretical, hitrate_theoretical, 'LineWidth', 1, 'Color', 'green')
%     plot(all_thresholds, ones(size(all_thresholds))/2, 'pg') % thresholds
%     hold off
%     % xlabel('SPL_{TH}')
%     xlabel('Signal ABR amplitude (\muV)', 'FontSize', 32)
%     ylabel('Probability of hit', 'FontSize', 32)
%     % title(['Psychometric curve for signal ABR detection (Vd=', num2str(Vd), ' {\muV}^2, FP=', num2str(FPR), ')'])
%     title(['Psychometric curve (Vd=', num2str(Vd), ' {\muV}^2)'], 'FontSize', 32)
%     box off % remove ticks on top and right borders of plot
%     % axis tight % makes edges of data flush with left and right borders of plot
%     legend(['New method (m=', num2str(m), ' single traces)'],'Old method (1 averaged trace)', ...
%         ['Theoretical (M=', num2str(m_theor), ' single traces)'], 'location','northeast', 'FontSize', 16)
%     legend boxoff % remove box around legend

    %%Plot psychometric-like curve
    h=figure; errorbar(A, hit_rates_new, error_bar_exp, 'LineWidth', 1, 'Color', 'red')
    hold on
    errorbar(A, hit_rates_old, error_bar_old, 'LineWidth', 1, 'Color', 'blue')
    errorbar(A_theoretical, hitrate_theoretical, error_bar_theor, 'LineWidth', 1, 'Color', 'green')
    plot(all_thresholds, ones(size(all_thresholds))/2, 'pg') % thresholds
    hold off
    % xlabel('SPL_{TH}')
    xlabel('Signal ABR amplitude (\muV)', 'FontSize', 32)
    ylabel('Probability of hit', 'FontSize', 32)
    % title(['Psychometric curve for signal ABR detection (Vd=', num2str(Vd), ' {\muV}^2, FP=', num2str(FPR), ')'])
    title(['Psychometric curve (Vd=', num2str(Vd), ' {\muV}^2)'], 'FontSize', 32)
    box off % remove ticks on top and right borders of plot
    % axis tight % makes edges of data flush with left and right borders of plot
    legend(['New method (m=', num2str(m), ' single traces)'],'Old method (1 averaged trace)', ...
        ['Theoretical (M=', num2str(m_theor), ' single traces)'], 'location','northeast', 'FontSize', 16)
    legend boxoff % remove box around legend

    % Save figure
    filename = ['psychometric_m', num2str(m), '_M', num2str(m_theor), '_Vd', num2str(Vd), '.fig'];
    saveas(h, filename)
    % close all

    err_Ath_old_cache(count) = err_Ath_old;
    err_Ath_new_cache(count) = err_Ath_new;
    Vd_cache(count) = Vd;
    Vd_old_cache(count) = Vd_old;
    V_0_cache(count) = V_0;
end

% 
% disp(SNR)
% disp(Fsp)
% disp(TH_old)
% disp(TH_new)