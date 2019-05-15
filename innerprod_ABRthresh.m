%% innerprod_ABRthresh.m
% 
% Calculate signal at each dB SPL level, using the inner product of single
% traces with (normalized) averaged ABR trace at highest dB level "signal
% basis". K-S or sine-rank test of inner product distributions with that at
% 0 dB SPL to calc p-value of signal at each dB level. 
% Method suggested by Stephen Bates
%
% Run this script after "import_ABRcsv_folder.m", which loads Noor's
% single-trace ABR data.
%
% Assesses: <signal(amp)> - averaged ABR voltage across all single-traces at the
% same time-point in each single trace.
%
% Dependencies: none
% Last edit: 5/13/2019
%
% Author: George Liu

SAMPLES = size(X_csv{1}, 1);
SAMPLING_RATE = 200000; % Hz
dt = 1/SAMPLING_RATE * 1000; % ms
T_total = (SAMPLES-1)*dt; % ms
x = 0:dt:T_total;

%% Plot averaved ABR waveforms -- the current standard for visually determining ABR threshold
% Also plot single-trace inner products

% Calculate signal basis vector from normalized max averaged ABR trace
max_averagedABR_trace = mean(X_csv{end}, 2); % SAMPLES x 1 vector
magn = dot(max_averagedABR_trace, max_averagedABR_trace);
signal_basis = max_averagedABR_trace / sqrt(magn); % SAMPLES x 1 vector

% Compute distribution of inner products (single traces) at each dB level
dist_innerprod = cell(A_length, 1);
for i = 1:A_length
    thisX = X_csv{i}; % SAMPLES x m_traces matrix
    % 5-12-19: Uncomment next 2 lines to normalize each single trace by its
    % norm squared, and change magn above to sqrt of its value
    thisX_norms = dot(thisX, thisX, 1); % 1 x m_traces matrix
    thisX = thisX./sqrt(thisX_norms); % SAMPLES x m_traces matrix
    dist_innerprod{i} = thisX' * signal_basis; % m_traces x 1 vector
end

% Plot averaged ABR trace and distribution of single-trace signal
% components per dB SPL level
count2 = 0;
figure('DefaultAxesFontSize', 16)
AxesHandles_3 = zeros(A_length, 1);
axesHandle = zeros(A_length, 1);
for i = 1:A_length
    % Plot averaged ABR trace in first column
    count2 = count2 + 1;
    AxesHandles_3(i) = subplot(A_length, 2, count2);  
    y = mean(X_csv{i}, 2);
    plot(x, y*10^6)
    title(['Average ABR, Input A=', num2str(A_csv(i)), ' dB SPL'])
    xlabel('Time (ms)')
    ylabel('\muV')
    
    % Plot histogram of single-trace inner products
    count2 = count2 + 1;
    axesHandle(i) = subplot(A_length, 2, count2); 
    histogram(dist_innerprod{i}, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'edgecolor', 'none')
    xlabel('Inner product (a.u.)')
    ylabel('Prob')
    title(['Input A=', num2str(A_csv(i)), ' dB SPL'])

    % plot mean as blue line
    hold on
    signal_mean = mean(dist_innerprod{i});
    line([signal_mean signal_mean], ylim, 'LineWidth', 1, 'Color', 'b'); % vertical line for cutoff
    hold off
end
same_yaxes(AxesHandles_3)
% Same x and y axes for all histograms
same_yaxes(axesHandle)
same_xaxes(axesHandle)

%% Calculate p-value for each dB SPL level, using Wilcoxin sign-rank test or K-S test
% compared with 0 dB SPL inner product distribution
zero_dB_innerdist = dist_innerprod{1};
p_val = zeros(A_length, 1);
signrank_stat = zeros(A_length, 1);
p_val_KS = zeros(A_length, 1);
ks_score = zeros(A_length,1);
for i = 1:A_length
    [p_val(i), ~, stats] = signrank(dist_innerprod{i}, zero_dB_innerdist); % Wilcoxin sign-rank test, two-sided 
    signrank_stat(i) = stats.signedrank;
    [~, p_val_KS(i), ks_score(i)] = kstest2(dist_innerprod{i}, zero_dB_innerdist); % Two-sample Kolmogorov-Smirnov test (unequal cdf's)
end

% Plot p-value vs dB SPL level: Wilcoxin sign rank test
figure('DefaultAxesFontSize', 20)
plot(A_csv, p_val, '-o')
title(['Wilcoxin sign rank for inner products vs 0 dB SPL'])
xlabel('Amplitude (dB SPL)')
ylabel('p-value')  
% CUSTOM add coordinate of threshold point
THRESH_INDEX = 5; % index of A_csv value that yields p_val < 0.05; from looking at plot
x_text = A_csv(THRESH_INDEX);
y_text = p_val(THRESH_INDEX);
text(x_text - 8, y_text, ['(' num2str(x_text) ', ' num2str(y_text, 2) ')'])

% Plot Wilcoxign sign rank z-stat vs dB SPL level: Wilcoxin sign rank test
figure('DefaultAxesFontSize', 20)
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
figure('DefaultAxesFontSize', 20)
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
figure('DefaultAxesFontSize', 20)
plot(A_csv, ks_score, '-o')
title(['Kolmogorov-Smirnov test for inner products vs 0 dB SPL'])
xlabel('Amplitude (dB SPL)')
ylabel('K-S score')  
% CUSTOM add coordinate of threshold point
THRESH_INDEX = 5; % index of A_csv value that yields p_val < 0.05; from looking at plot
x_text = A_csv(THRESH_INDEX);
y_text = ks_score(THRESH_INDEX);
text(x_text - 8, y_text, ['(' num2str(x_text) ', ' num2str(y_text, 2) ')'])

%% Compute cross-correlation instead of inner products, to add shifts for latency

% Calculate signal basis vector from normalized max averaged ABR trace
max_averagedABR_trace = mean(X_csv{end}, 2); % SAMPLES x 1 vector
magn = dot(max_averagedABR_trace, max_averagedABR_trace);
signal_basis = max_averagedABR_trace / magn; % SAMPLES x 1 vector

% Compute distribution of cross correlation maxima (single traces) at each dB level
MAXLAG = 0.8/dt; % hard-code max lag of 1 ms, 5-14-19
dist_crosscor = cell(A_length, 1);
dist_lags = cell(A_length, 1);
for i = 1:A_length
    thisX = X_csv{i}; % SAMPLES x m_traces matrix
    m_traces = size(thisX, 2);
    max_xcorr = zeros(m_traces, 1);
    max_xcorr_lag = zeros(m_traces, 1);
    for j = 1:m_traces
        [r,lags] = xcorr(thisX(:, j), signal_basis, MAXLAG);
        [max_xcorr(j), ind] = max(r);
        max_xcorr_lag(j) = lags(ind);
    end
    dist_crosscor{i} = max_xcorr; % m_traces x 1 vector
    dist_lags{i} = max_xcorr_lag*dt; % m_traces x 1 vector, change units to ms
end

% Plot averaged ABR trace and distribution of single-trace signal
% components per dB SPL level
count2 = 0;
figure('DefaultAxesFontSize', 16)
AxesHandles_3 = zeros(A_length, 1);
axesHandle = zeros(A_length, 1);
axesHandle_2 = zeros(A_length, 1);
for i = 1:A_length
    % Plot averaged ABR trace in first column
    count2 = count2 + 1;
    AxesHandles_3(i) = subplot(A_length, 3, count2);  
    y = mean(X_csv{i}, 2);
    plot(x, y*10^6)
    title(['Average ABR, Input A=', num2str(A_csv(i)), ' dB SPL'])
    xlabel('Time (ms)')
    ylabel('\muV')
    
    % Plot histogram of single-trace cross-correlation maxima
    count2 = count2 + 1;
    axesHandle(i) = subplot(A_length, 3, count2); 
    histogram(dist_crosscor{i}, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'edgecolor', 'none')
    xlabel('Max xcorr (a.u.)')
    ylabel('Prob')
    title(['Input A=', num2str(A_csv(i)), ' dB SPL'])

    % plot mean as blue line
    hold on
    signal_mean = mean(dist_crosscor{i});
    line([signal_mean signal_mean], ylim, 'LineWidth', 1, 'Color', 'b'); % vertical line for cutoff
    hold off
    
    % Plot histogram of single-trace cross-correlation maxima lags
    count2 = count2 + 1;
    axesHandle_2(i) = subplot(A_length, 3, count2); 
    thislags_ms = dist_lags{i};  % ms units
    histogram(thislags_ms, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'edgecolor', 'none')
    xlabel('Lag (ms)')
    ylabel('Prob')
%     title(['Input A=', num2str(A_csv(i)), ' dB SPL'])

    % plot mean as blue line
    hold on
    signal_mean = mean(thislags_ms);
    line([signal_mean signal_mean], ylim, 'LineWidth', 1, 'Color', 'r'); % vertical line for cutoff
    hold off
end
same_yaxes(AxesHandles_3)
% Same x and y axes for all histograms
same_yaxes(axesHandle)
same_xaxes(axesHandle)
same_yaxes(axesHandle_2)
same_xaxes(axesHandle_2)

%% Calculate p-value for each dB SPL level, using Wilcoxin sign-rank test or K-S test
% compared with 0 dB SPL cross correlation distribution
zero_dB_innerdist = dist_crosscor{1};
p_val = zeros(A_length, 1);
signrank_stat = zeros(A_length, 1);
p_val_KS = zeros(A_length, 1);
ks_score = zeros(A_length,1);
for i = 1:A_length
    [p_val(i), ~, stats] = signrank(dist_crosscor{i}, zero_dB_innerdist); % Wilcoxin sign-rank test, two-sided 
    signrank_stat(i) = stats.signedrank;
    [~, p_val_KS(i), ks_score(i)] = kstest2(dist_crosscor{i}, zero_dB_innerdist); % Two-sample Kolmogorov-Smirnov test (unequal cdf's)
end

% Plot p-value vs dB SPL level: Wilcoxin sign rank test
figure('DefaultAxesFontSize', 20)
plot(A_csv, p_val, '-o')
title(['Wilcoxin sign rank for cross correlation maxima vs 0 dB SPL'])
xlabel('Amplitude (dB SPL)')
ylabel('p-value')  
% CUSTOM add coordinate of threshold point
THRESH_INDEX = 5; % index of A_csv value that yields p_val < 0.05; from looking at plot
x_text = A_csv(THRESH_INDEX);
y_text = p_val(THRESH_INDEX);
text(x_text - 8, y_text, ['(' num2str(x_text) ', ' num2str(y_text, 2) ')'])

% Plot Wilcoxign sign rank z-stat vs dB SPL level: Wilcoxin sign rank test
figure('DefaultAxesFontSize', 20)
plot(A_csv, signrank_stat, '-o')
title(['Wilcoxin sign rank for xcorr maxima vs 0 dB SPL'])
xlabel('Amplitude (dB SPL)')
ylabel('Sign rank test statistic')  
% CUSTOM add coordinate of threshold point
THRESH_INDEX = 5; % index of A_csv value that yields p_val < 0.05; from looking at plot
x_text = A_csv(THRESH_INDEX);
y_text = signrank_stat(THRESH_INDEX);
text(x_text - 8, y_text, ['(' num2str(x_text) ', ' num2str(y_text, 2) ')'])

% Plot p-value vs dB SPL level: K-S test
figure('DefaultAxesFontSize', 20)
plot(A_csv, p_val_KS, '-o')
title(['Kolmogorov-Smirnov test for xcorr maxima vs 0 dB SPL'])
xlabel('Amplitude (dB SPL)')
ylabel('p-value')  
% CUSTOM add coordinate of threshold point
THRESH_INDEX = 5; % index of A_csv value that yields p_val < 0.05; from looking at plot
x_text = A_csv(THRESH_INDEX);
y_text = p_val_KS(THRESH_INDEX);
text(x_text - 8, y_text, ['(' num2str(x_text) ', ' num2str(y_text, 2) ')'])

% Plot KS_score vs dB SPL level: K-S test
figure('DefaultAxesFontSize', 20)
plot(A_csv, ks_score, '-o')
title(['Kolmogorov-Smirnov test for xcorr maxima vs 0 dB SPL'])
xlabel('Amplitude (dB SPL)')
ylabel('K-S score')  
% CUSTOM add coordinate of threshold point
THRESH_INDEX = 5; % index of A_csv value that yields p_val < 0.05; from looking at plot
x_text = A_csv(THRESH_INDEX);
y_text = ks_score(THRESH_INDEX);
text(x_text - 8, y_text, ['(' num2str(x_text) ', ' num2str(y_text, 2) ')'])