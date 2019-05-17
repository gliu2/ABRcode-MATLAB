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

%% Inner products: Plot averaved ABR waveforms -- the current standard for visually determining ABR threshold
% Also plot single-trace inner products

% Calculate signal basis vector from normalized max averaged ABR trace
max_averagedABR_trace = mean(X_csv{end}, 2); % SAMPLES x 1 vector
magn = dot(max_averagedABR_trace, max_averagedABR_trace);
signal_basis = max_averagedABR_trace / magn; % SAMPLES x 1 vector
% signal_basis = max_averagedABR_trace / sqrt(magn); % SAMPLES x 1 vector

% Compute distribution of inner products (single traces) at each dB level
dist_innerprod = cell(A_length, 1);
for i = 1:A_length
    thisX = X_csv{i}; % SAMPLES x m_traces matrix
    % 5-12-19: Uncomment next 2 lines to normalize each single trace by its
    % norm squared, and change magn above to sqrt of its value
%     thisX_norms = dot(thisX, thisX, 1); % 1 x m_traces matrix
%     thisX = thisX./sqrt(thisX_norms); % SAMPLES x m_traces matrix
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

%% 5-15-19: Find latency lags of averaged ABR traces
% Calculate signal basis vector from normalized max averaged ABR trace
max_averagedABR_trace = mean(X_csv{end}, 2); % SAMPLES x 1 vector
magn = sqrt(dot(max_averagedABR_trace, max_averagedABR_trace));
signal_basis = max_averagedABR_trace/magn; % SAMPLES x 1 vector

% Find lags using cross-correlation maxima of averaged ABR traces
avg_lag = zeros(A_length, 1);
for i = 1:A_length
    % Plot averaged ABR trace in first column
    y = mean(X_csv{i}, 2);
    [r,lags] = xcorr(y, signal_basis);
    [~, ind] = max(r);
    avg_lag(i) = lags(ind)*dt;
end

figure('DefaultAxesFontSize', 20)
plot(A_csv, avg_lag)
xlabel('Stimulus level (dB SPL)')
ylabel('Lag (ms)')
title('Average ABR cross-correlation maxima lags')

%% 5-15-19: Find lags using peak detection method (PTDetect.m) on relative scale
% averaged ABRs
PEAK_THRESH = 0.5;
yrel_cache = cell(A_length, 1);
peak_cache = cell(A_length, 1);
trough_cache = cell(A_length, 1);

figure
AxesHandles_1 = zeros(A_length, 1);
for i = 1:A_length
    y = mean(X_csv{i}, 2);
    y_rel = y/max(abs(y)); % y vector is scaled to max 1
    
    AxesHandles_1(i) = subplot(ceil(A_length/2),2,i);
    plot(x, y_rel)
    title(['Input A=', num2str(A_csv(i)), ' dB SPL'])
    xlabel('Time (ms)')
    ylabel('Relative voltage (a.u.)')
    
    % Plot peaks and troughs
    hold on
    [P, T] = PTDetect(y_rel, PEAK_THRESH);
    scatter(x(P), y_rel(P), 'b') % blue peaks
    scatter(x(T), y_rel(T), 'r') % red troughs
%     hold off
    
    % cache variables
    yrel_cache{i} = y_rel;
    peak_cache{i} = P;
    trough_cache{i} = T;
end
same_yaxes(AxesHandles_1)

%TODO: trace peak of big wave, using nearest neighbor across dB to measure
%latency. Then use latencies from both methods to (1) correct and inner
%product, (2) cross-correlate with restricted lag

% Nearest peak to measure latency
maxpeaks_ind = zeros(A_length, 1);
[~, ind] = max(yrel_cache{end}); % find max peak at highest level ABR ave
disp('Is largest peak detected:')
disp(any(P==ind))
maxpeaks_ind(end) = ind;
for i = A_length-1:-1:1
    lower_peaks = peak_cache{i};
    [~, i_peak] = min(abs(lower_peaks - maxpeaks_ind(i+1)));
    maxpeaks_ind(i) = lower_peaks(i_peak);
end

% Mark max peak in previous plots
for i=1:A_length
    scatter(x(maxpeaks_ind(i)),  yrel_cache{i}(maxpeaks_ind(i)), 'g', 'Parent', AxesHandles_1(i))
end

%% Plot lag versus level calculated with peak-measured shifts
avg_lag2 = (x(maxpeaks_ind) - x(maxpeaks_ind(end)));

figure('DefaultAxesFontSize', 20)
plot(A_csv, avg_lag2, '-ob')
xlabel('Stimulus level (dB SPL)')
ylabel('Lag (ms)')
title('Average ABR max peak lags')

%% 5-16-19 Recalculate lags of peaks using only peak-positions of 50 db SPL and
% greater (less may be unreliable), linear regression
midA_ind = round(A_length/2);
XX = A_csv - A_csv(end); % A_length x 1 vector, shifted to zero-intercept (zero lag) at max amplitude
b1 = XX(midA_ind:end)\avg_lag2(midA_ind:end)'; % b1 = x\y; Find the linear regression relation y=?1x between amplitudes x (db SPL) and lags y (ms), uses least squares
% mdl = fitlm(A,b,'linear','Intercept',false)%'RobustOpts','off' is the default
lag0 = -A_csv(end)*b1;
avg_lag3 = lag0 + b1*A_csv;
% Plot lag versus level calculated with peak-measured shifts, linear
% regression method
% figure('DefaultAxesFontSize', 20)
hold on
plot(A_csv, avg_lag3, '--g')
% xlabel('Stimulus level (dB SPL)')
% ylabel('Lag (ms)')
hold off
title(['Peak lags with fit for ', num2str(A_csv(midA_ind)), ' dB peaks and above'])

%% 5-16-19: Inner products of single-trace inner products with latency-corrected average ABR max peak
% latencies calculated using peak-lag method above
% THEN RUN SIGN RANK AND K-S TESTS USING CELL WAY ABOVE

% Calculate signal basis vector from normalized max averaged ABR trace
max_averagedABR_trace = mean(X_csv{end}, 2); % SAMPLES x 1 vector
magn = dot(max_averagedABR_trace, max_averagedABR_trace);
signal_basis = max_averagedABR_trace / magn; % SAMPLES x 1 vector

% Compute distribution of inner products (single traces) at each dB level
dist_innerprod = cell(A_length, 1);
for i = 1:A_length
    thisX = X_csv{i}; % SAMPLES x m_traces matrix
    shifted_signalbasis = zeros(size(signal_basis));
    ind_shift = round(avg_lag2(i)/dt); % lags from peak shifts at all dB level
%     ind_shift = round(avg_lag3(i)/dt); % use lags from fitting max peaks lags at 45 dB SPL and higher, ignoring lower db peaks which are less reliable 
    shifted_signalbasis(1 + ind_shift:end) = signal_basis(1:end - ind_shift);
    dist_innerprod{i} = thisX' * shifted_signalbasis; % m_traces x 1 vector
end

% Plot averaged ABR trace and distribution of single-trace signal
% components per dB SPL level
count2 = 0;
figure('DefaultAxesFontSize', 16)
AxesHandles_3 = zeros(A_length, 1);
AxesHandles_4 = zeros(A_length, 1);
axesHandle = zeros(A_length, 1);
for i = 1:A_length
    % Plot averaged ABR trace in first column - absolute scale
    count2 = count2 + 1;
    AxesHandles_3(i) = subplot(A_length, 3, count2);  
    y = mean(X_csv{i}, 2);
    plot(x, y*10^6)
    title(['Average ABR, Input A=', num2str(A_csv(i)), ' dB SPL'])
    xlabel('Time (ms)')
    ylabel('\muV')
    
    % Plot averaged ABR trace in second column - relative scale
    count2 = count2 + 1;
    AxesHandles_4(i) = subplot(A_length, 3, count2);  
%     y = mean(X_csv{i}, 2);
%     plot(x, y*10^6)
    plot(x, yrel_cache{i})
    title(['Average ABR, Input A=', num2str(A_csv(i)), ' dB SPL'])
    xlabel('Time (ms)')
    ylabel('Rel volt')
    
    % Plot peaks and troughs as overlay on relative average ABR
    hold on
    scatter(x(peak_cache{i}), yrel_cache{i}(peak_cache{i}), 'b') % blue peaks
    scatter(x(trough_cache{i}), yrel_cache{i}(trough_cache{i}), 'r') % red troughs
    scatter(x(maxpeaks_ind(i)),  yrel_cache{i}(maxpeaks_ind(i)), 'g')
    hold off
    
    % Plot histogram of single-trace inner products
    count2 = count2 + 1;
    axesHandle(i) = subplot(A_length, 3, count2); 
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
% Same y axes for average ABR plots
same_yaxes(AxesHandles_3)
same_yaxes(AxesHandles_4)
% Same x and y axes for all histograms
same_yaxes(axesHandle)
same_xaxes(axesHandle)


%% 5-17-19: Chunk inner product: histograms of Inner products like above, but use trough-trough chunk around max peak to do lag-inner product
maxdb_troughrelmaxpeak = trough_cache{end} - maxpeaks_ind(end);
chunk_start = max(maxdb_troughrelmaxpeak(maxdb_troughrelmaxpeak<0)) + maxpeaks_ind(end);
chunk_end = min(maxdb_troughrelmaxpeak(maxdb_troughrelmaxpeak>0)) + maxpeaks_ind(end);
% [~, ind] = max(maxdb_troughrelmaxpeak(maxdb_troughrelmaxpeak<0));
% chunk_start = trough_cache{end}(ind);
% [~, ind] = min(maxdb_troughrelmaxpeak(maxdb_troughrelmaxpeak>0));
% chunk_end = trough_cache{end}(ind);

% Calculate signal basis vector from normalized max averaged ABR trace
max_averagedABR_trace = mean(X_csv{end}, 2); % SAMPLES x 1 vector
% Additional step: Convert signal basis to chunk
max_averagedABR_trace(1:chunk_start-1) = 0; 
max_averagedABR_trace(chunk_end+1:end) = 0; 
% Normalize
magn = dot(max_averagedABR_trace, max_averagedABR_trace);
signal_basis = max_averagedABR_trace / magn; % SAMPLES x 1 vector


% Compute distribution of inner products (single traces) at each dB level
dist_innerprod = cell(A_length, 1);
for i = 1:A_length
    thisX = X_csv{i}; % SAMPLES x m_traces matrix
    shifted_signalbasis = zeros(size(signal_basis));
    ind_shift = round(avg_lag2(i)/dt); % lags from peak shifts at all dB level
%     ind_shift = round(avg_lag3(i)/dt); % use lags from fitting max peaks lags at 45 dB SPL and higher, ignoring lower db peaks which are less reliable 
    shifted_signalbasis(1 + ind_shift:end) = signal_basis(1:end - ind_shift);
    dist_innerprod{i} = thisX' * shifted_signalbasis; % m_traces x 1 vector
end

% Plot averaged ABR trace and distribution of single-trace signal
% components per dB SPL level
count2 = 0;
figure('DefaultAxesFontSize', 16)
AxesHandles_3 = zeros(A_length, 1);
AxesHandles_4 = zeros(A_length, 1);
axesHandle = zeros(A_length, 1);
for i = 1:A_length
    % Plot averaged ABR trace in first column - absolute scale
    count2 = count2 + 1;
    AxesHandles_3(i) = subplot(A_length, 3, count2);  
    y = mean(X_csv{i}, 2);
    plot(x, y*10^6)
    title(['Average ABR, Input A=', num2str(A_csv(i)), ' dB SPL'])
    xlabel('Time (ms)')
    ylabel('\muV')
    
    % Plot averaged ABR trace in second column - relative scale
    count2 = count2 + 1;
    AxesHandles_4(i) = subplot(A_length, 3, count2);  
%     y = mean(X_csv{i}, 2);
%     plot(x, y*10^6)
    plot(x, yrel_cache{i})
    title(['Average ABR, Input A=', num2str(A_csv(i)), ' dB SPL'])
    xlabel('Time (ms)')
    ylabel('Rel volt')
    
    % Plot peaks and troughs as overlay on relative average ABR
    hold on
    scatter(x(peak_cache{i}), yrel_cache{i}(peak_cache{i}), 'b') % blue peaks
    scatter(x(trough_cache{i}), yrel_cache{i}(trough_cache{i}), 'r') % red troughs
    scatter(x(maxpeaks_ind(i)),  yrel_cache{i}(maxpeaks_ind(i)), 'g')
    hold off
    
    % Plot vertical lines denoting max average ABR "chunk" bounds at
    % troughs
    if i==A_length
        hold on
        line([x(chunk_start) x(chunk_start)], ylim, 'LineWidth', 1, 'Color', 'k'); % vertical line for chunk start
        line([x(chunk_end) x(chunk_end)], ylim, 'LineWidth', 1, 'Color', 'k'); % vertical line for chunk end
        hold off
    end
    
    % Plot histogram of single-trace inner products
    count2 = count2 + 1;
    axesHandle(i) = subplot(A_length, 3, count2); 
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
% Same y axes for average ABR plots
same_yaxes(AxesHandles_3)
same_yaxes(AxesHandles_4)
% Same x and y axes for all histograms
same_yaxes(axesHandle)
same_xaxes(axesHandle)

%% Compute cross-correlation instead of inner products, to add shifts for latency

% Calculate signal basis vector from normalized max averaged ABR trace
max_averagedABR_trace = mean(X_csv{end}, 2); % SAMPLES x 1 vector
magn = dot(max_averagedABR_trace, max_averagedABR_trace); % normalize by average ABR norm squared
% magn = sqrt(dot(max_averagedABR_trace, max_averagedABR_trace)); % normalize by single trace norm and average ABR norm
signal_basis = max_averagedABR_trace / magn; % SAMPLES x 1 vector

% Compute distribution of cross correlation maxima (single traces) at each dB level
MAXLAG = 0.8/dt; % hard-code max lag of 1 ms, 5-14-19
dist_crosscor = cell(A_length, 1);
dist_lags = cell(A_length, 1);
for i = 1:A_length
    thisX = X_csv{i}; % SAMPLES x m_traces matrix
    thisX_norms = dot(thisX, thisX, 1); % 1 x m_traces matrix
    m_traces = size(thisX, 2);
    max_xcorr = zeros(m_traces, 1);
    max_xcorr_lag = zeros(m_traces, 1);
    for j = 1:m_traces
        [r,lags] = xcorr(thisX(:, j), signal_basis, MAXLAG); % normalize by average ABR norm squared
%         [r,lags] = xcorr(thisX(:, j)/thisX_norms(j), signal_basis, MAXLAG); % normalize by single trace norm and average ABR norm
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