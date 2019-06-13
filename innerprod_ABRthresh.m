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
% Dependencies: same_yaxes.m, same_xaxes.m, PTDetect.m,
%               analyze_innerprod_ABR.m, innerprod2pval.m, visualize_innerprod_stats.m
%
% Last edit: 6-8/2019 - modularized calculating innerprods, pvals,
%                        visualize stats plots
%
% Author: George Liu

SAMPLES = size(X_csv{1}, 1);
SAMPLING_RATE = 200000; % Hz
dt = 1/SAMPLING_RATE * 1000; % ms
T_total = (SAMPLES-1)*dt; % ms
x = 0:dt:T_total;

%% Inner products: Plot averaved ABR waveforms -- the current standard for visually determining ABR threshold
% Also plot single-trace inner products

dist_innerprod = analyze_innerprod_ABR(X_csv, zeros(A_length, 1));
% % Calculate signal basis vector from normalized max averaged ABR trace
% max_averagedABR_trace = mean(X_csv{end}, 2); % SAMPLES x 1 vector
% magn = dot(max_averagedABR_trace, max_averagedABR_trace);
% signal_basis = max_averagedABR_trace / magn; % SAMPLES x 1 vector
% % signal_basis = max_averagedABR_trace / sqrt(magn); % SAMPLES x 1 vector
% 
% % Compute distribution of inner products (single traces) at each dB level
% dist_innerprod = cell(A_length, 1);
% for i = 1:A_length
%     thisX = X_csv{i}; % SAMPLES x m_traces matrix
%     % 5-12-19: Uncomment next 2 lines to normalize each single trace by its
%     % norm squared, and change magn above to sqrt of its value
% %     thisX_norms = dot(thisX, thisX, 1); % 1 x m_traces matrix
% %     thisX = thisX./sqrt(thisX_norms); % SAMPLES x m_traces matrix
%     dist_innerprod{i} = thisX' * signal_basis; % m_traces x 1 vector
% end

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

[p_val, signrank_stat, p_val_KS, ks_score] = innerprod2pval(dist_innerprod); % 5-26-2019 implementation
visualize_innerprod_stats(p_val, signrank_stat, p_val_KS, ks_score, A_csv)

%% 5-25-19: Find latency lags of averaged ABR traces using cross-covariance with coeff normalization

% Calculate signal basis vector from normalized max averaged ABR trace
max_averagedABR_trace = mean(X_csv{end}, 2); % SAMPLES x 1 vector
magn = sqrt(dot(max_averagedABR_trace, max_averagedABR_trace));
signal_basis = max_averagedABR_trace/magn; % SAMPLES x 1 vector

% Find lags using cross-correlation maxima of averaged ABR traces
avg_lag = zeros(A_length, 1);
r_cache = cell(A_length, 1);
lags_cache = cell(A_length, 1);
maxr_cache = zeros(A_length, 1);
for i = 1:A_length
    % Plot averaged ABR trace in first column
    y = mean(X_csv{i}, 2);
    [r,lags] = xcov(y, signal_basis, 'coeff');
    [max_r, ind] = max(r);
    avg_lag(i) = lags(ind)*dt;
    
    % cache vars
    r_cache{i} = r;
    lags_cache{i} = lags;
    maxr_cache(i) = max_r;
end

% Plot xcov lags versus dB level
hh=figure('DefaultAxesFontSize', 20);
plot(A_csv, avg_lag, 'o-b');
xlabel('Stimulus level (dB SPL)')
ylabel('Lag (ms)')
title('Average ABR cross-covariance maxima lags')

% 5-25-19: Plot ave ABR cross-correlations (for finding lags) with signal basis to identify good peak
%locations
figure
AxesHandles_1 = zeros(A_length, 1);
for i = 1:A_length
    AxesHandles_1(i) = subplot(ceil(A_length/2),2,i);
    plot(lags_cache{i}, r_cache{i})
    title(['Input A=', num2str(A_csv(i)), ' dB SPL'])
    xlabel('Time (ms)')
    ylabel('xcov (a.u.)')
end
same_yaxes(AxesHandles_1)
    
% Extrapolate starting from lowest dB where max cross-covariance is > 0.50
% Find lowest dB
XCOV_THRESH_LAG = 0.5; % for coeff normalized cross-covariance
[row, col] = find(maxr_cache > XCOV_THRESH_LAG);
ind_uselag = min(row);

% Extrapolate
avg_lag_extrap = avg_lag(ind_uselag:end);
avg_lag_extrap = interp1(avg_lag_extrap, -(ind_uselag-2):1, 'linear', 'extrap')';
figure(hh)
hold on
plot(A_csv(1:length(avg_lag_extrap)), avg_lag_extrap, 'o--g')
title(['Cross-covariance lags, max(xcov)>', num2str(XCOV_THRESH_LAG)])
hold off

%TODO 5-27-19: Ensure that extrapolated lags have ceiling of max possible lag

% Plot maximum cross-covariance versus dB level
figure('DefaultAxesFontSize', 20)
plot(A_csv, maxr_cache, 'o-b')
hold on, yline(XCOV_THRESH_LAG); hold off % add cutoff line to max xcov vs dB plot
xlabel('Stimulus level (dB SPL)')
ylabel('Max xcov (a.u.)')
title(['Max cross-covariance, cutoff=', num2str(XCOV_THRESH_LAG)])

% Save final extrapolation-correct xcov lags
avg_lag_xcovExtrap = avg_lag;
avg_lag_xcovExtrap(1:ind_uselag-1) = avg_lag_extrap(1:ind_uselag-1);


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
% 5-25-19 added to modularize lag calculation
% 5-16-19: latencies calculated using peak-lag method above
% THEN RUN SIGN RANK AND K-S TESTS USING CELL WAY ABOVE
this_lag = avg_lag_xcovExtrap; % 5-25-19 added to modularize lag calculation; this_lag units of ms

dist_innerprod = analyze_innerprod_ABR(X_csv, round(this_lag/dt));
% % Calculate signal basis vector from normalized max averaged ABR trace
% max_averagedABR_trace = mean(X_csv{end}, 2); % SAMPLES x 1 vector
% magn = dot(max_averagedABR_trace, max_averagedABR_trace);
% % signal_basis = max_averagedABR_trace / magn; % SAMPLES x 1 vector
% signal_basis = max_averagedABR_trace / sqrt(magn); % SAMPLES x 1 vector - 5-25-19: manual coeff normalization
% 
% % Compute distribution of inner products (single traces) at each dB level
% dist_innerprod = cell(A_length, 1);
% for i = 1:A_length
%     thisX = X_csv{i}; % SAMPLES x m_traces matrix
%     thisX = thisX./sqrt(var(thisX)); % 5-25-19: normalize single trace by standard deviation, column by column
%     shifted_signalbasis = zeros(size(signal_basis));
%     ind_shift = round(this_lag(i)/dt); % lags from peak shifts at all dB level
% %     ind_shift = round(avg_lag2(i)/dt); % lags from peak shifts at all dB level
% %     ind_shift = round(avg_lag3(i)/dt); % use lags from fitting max peaks lags at 45 dB SPL and higher, ignoring lower db peaks which are less reliable 
%     shifted_signalbasis(1 + ind_shift:end) = signal_basis(1:end - ind_shift);
%     dist_innerprod{i} = thisX' * shifted_signalbasis; % m_traces x 1 vector
% end

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
    
%     % 5-25-19: PEAK LAGS: Plot peaks and troughs as overlay on relative average ABR
%     hold on
%     scatter(x(peak_cache{i}), yrel_cache{i}(peak_cache{i}), 'b') % blue peaks
%     scatter(x(trough_cache{i}), yrel_cache{i}(trough_cache{i}), 'r') % red troughs
%     scatter(x(maxpeaks_ind(i)),  yrel_cache{i}(maxpeaks_ind(i)), 'g')
%     hold off

    % 5-25-19: XCOV LAGS: plot max peak lag as black line
    hold on
    peaklag_i = x(maxpeaks_ind(A_length))+this_lag(i);
    line([peaklag_i peaklag_i], ylim, 'LineWidth', 1, 'Color', 'k'); % vertical line for cutoff
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



% %% 5-27-19: Check if p-val estimates for lag-inner product method are robust to 0 dB traces random lag
% m0_traces = size(X_csv{1}, 2);
% lag0 = -(m0_traces-1):(m0_traces-1); % all possible lags for lag0; index, not ms
% orig_lags = this_lag;
% pval_wsr_all = zeros(A_length, length(lag0));
% stat_wsr_all = zeros(A_length, length(lag0));
% pval_ks_all = zeros(A_length, length(lag0));
% stat_ks_all = zeros(A_length, length(lag0));
% 
% for i = 1:length(lag0)
%     disp(['Working on ', num2str(i), ' out of ', num2str(length(lag0)), '....'])
%     lag_i = round(orig_lags/dt);
%     lag_i(1) = lag0(i); % test 0 dB lag
%     dist_innerprod_i = analyze_innerprod_ABR(X_csv, lag_i);
%     [pval_wsr_all(:, i), stat_wsr_all(:, i), pval_ks_all(:, i), stat_ks_all(:, i)] = innerprod2pval(dist_innerprod_i); 
% end
%     
% % Plot p-vals
% figure('DefaultAxesFontSize', 20)
% figure, plot(A_csv, pval_wsr_all)
% xlabel('Stimulus level (dB)')
% ylabel('P-val')
% title('Wilcoxin sign rank test, 0 dB lag 0 to +/-2.56 ms')

%% 5-17-19: Chunk inner product: histograms of Inner products like above, but use trough-trough chunk around max peak to do lag-inner product
% 5-22-19: Make chunk from t=0 to trough after peak
maxdb_troughrelmaxpeak = trough_cache{end} - maxpeaks_ind(end);
t_endnoise = find(x==2);
chunk_start = t_endnoise;
% chunk_start = max(maxdb_troughrelmaxpeak(maxdb_troughrelmaxpeak<0)) + maxpeaks_ind(end);  % Uncomment to make chunk limited to between troughs around peak
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

[p_val, signrank_stat, p_val_KS, ks_score] = innerprod2pval(dist_crosscor); % 5-26-2019 implementation
visualize_innerprod_stats(p_val, signrank_stat, p_val_KS, ks_score, A_csv) %Change figure titles to "xcorr maxima"

%% 5-23-19: Find lags of 3 chunks- (1) 2ms (stimulus onset) to 1st trough, (2) max peak, (3) 2nd-3rd trough
% uses peak detection method (PTDetect.m) on relative scale averaged ABRs
% Compare with 5-15-19 code
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

%% 5-23-19: Trace troughs to left and right of max peak to identify other chunk
% latencies
lefttroughs_ind = zeros(A_length, 1);
righttroughs_ind = zeros(A_length, 1);

% Get max dB level troughs locations
trough_dist = trough_cache{end}-maxpeaks_ind(end);
lefttroughs_ind(end) = trough_cache{end}(find(trough_dist==max(trough_dist(trough_dist<0))));
righttroughs_ind(end) = trough_cache{end}(find(trough_dist==min(trough_dist(trough_dist>0))));

% Trace lower dB level troughs
for i = A_length-1:-1:1
    lower_troughs = trough_cache{i};
    
    % left trough
    [~, i_trough] = min(abs(lower_troughs - lefttroughs_ind(i+1)));
    lefttroughs_ind(i) = lower_troughs(i_trough);
    
    % right trough
    [~, i_trough] = min(abs(lower_troughs - righttroughs_ind(i+1)));
    righttroughs_ind(i) = lower_troughs(i_trough);
end

% Mark troughs in previous plots
for i=1:A_length
    scatter(x(lefttroughs_ind(i)),  yrel_cache{i}(lefttroughs_ind(i)), 'c', 'Parent', AxesHandles_1(i))
    scatter(x(righttroughs_ind(i)),  yrel_cache{i}(righttroughs_ind(i)), 'm', 'Parent', AxesHandles_1(i))
end

%% 5-23-19: Plot lag versus level calculated with peak-measured shifts
maxpeak_lags = (x(maxpeaks_ind) - x(maxpeaks_ind(end))); % max peak lag (ms)
lefttrough_lags = (x(lefttroughs_ind) - x(lefttroughs_ind(end)));
righttrough_lags = (x(righttroughs_ind) - x(righttroughs_ind(end)));

figure('DefaultAxesFontSize', 20)
plot(A_csv, maxpeak_lags, '-ob')
hold on
plot(A_csv, lefttrough_lags, '-oc')
plot(A_csv, righttrough_lags, '-om')
hold off
xlabel('Stimulus level (dB SPL)')
ylabel('Lag (ms)')
title('Average ABR max peak lags')
legend({'Max peak', 'Trough 1', 'Trough 2'})

%% 5-25-19: Experiment if set all 3 chunk lags to that calculated by xcov max > 0.45, linear fit
maxpeak_lags = avg_lag_xcovExtrap; % max peak lag (ms)
lefttrough_lags = avg_lag_xcovExtrap;
righttrough_lags = avg_lag_xcovExtrap;

%% 5-23-19: Multiple plots of Chunk inner product: histograms of Inner products like above, but use trough-trough chunk around max peak to do lag-inner product
% Make chunks: (1) 2ms (stimulus onset) to 1st trough, (2) max peak, (3) 2nd-3rd trough

% Find 3rd trough (2nd after max peak)
troughs_dist2righttrough = trough_cache{end} - righttroughs_ind(end);
troughs2right = troughs_dist2righttrough(troughs_dist2righttrough > 0);
if ~any(troughs2right)
    disp('Warning: No 3rd trough found in ave ABR max dB level!')
end
trough3_ind = min(troughs2right) + righttroughs_ind(end);

% Get chunk start and end indices
t_endnoise = find(x==2);
chunk1_start = t_endnoise; % t=2ms (stimulus onset) to 1st trough
chunk1_end = lefttroughs_ind(end);
chunk2_start = chunk1_end; % 1st-2nd trough
chunk2_end = righttroughs_ind(end);
chunk3_start = chunk2_end; % 2nd-3rd trough
chunk3_end = trough3_ind;

% Calculate signal basis vector from normalized max averaged ABR trace
max_averagedABR_trace = mean(X_csv{end}, 2); % SAMPLES x 1 vector
% Additional step: Convert signal basis to chunks
% chunk 1: t=2ms (stimulus onset) to 1st trough
chunk1 = max_averagedABR_trace;
chunk1(1:chunk1_start-1) = 0;
chunk1(chunk1_end+1:end) = 0; 

% chunk 2: 1st-2nd trough
chunk2 = max_averagedABR_trace;
chunk2(1:chunk2_start-1) = 0;
chunk2(chunk2_end+1:end) = 0; 

% chunk 3: 2nd-3rd trough
chunk3 = max_averagedABR_trace;
chunk3(1:chunk3_start-1) = 0;
chunk3(chunk3_end+1:end) = 0; 

% Compute distribution of inner products (single traces) at each dB level. % m_traces x 1 vector
dist_innerprod1 = analyze_innerprod_ABR(X_csv, round(lefttrough_lags/dt), chunk1); % 1st trough lag; lags from peak shifts at all dB level
dist_innerprod2 = analyze_innerprod_ABR(X_csv, round(maxpeak_lags/dt), chunk2); % Max peak lags; lags from peak shifts at all dB level
dist_innerprod3 = analyze_innerprod_ABR(X_csv, round(righttrough_lags/dt), chunk3); % 2nd trough lag; lags from peak shifts at all dB level

% Plot averaged ABR trace and distribution of single-trace signal
% components per dB SPL level
count2 = 0;
figure('DefaultAxesFontSize', 16)
AxesHandles_4 = zeros(A_length, 1);
axesHandle1 = zeros(A_length, 1);
axesHandle2 = zeros(A_length, 1);
axesHandle3 = zeros(A_length, 1);
for i = 1:A_length
    % Plot averaged ABR trace in second column - relative scale
    count2 = count2 + 1;
    AxesHandles_4(i) = subplot(A_length, 4, count2);  
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
    scatter(x(maxpeaks_ind(i)),  yrel_cache{i}(maxpeaks_ind(i)), 'g') % max peak in green
    scatter(x(lefttroughs_ind(i)),  yrel_cache{i}(lefttroughs_ind(i)), 'c') % trough to left of max peak in cyan
    scatter(x(righttroughs_ind(i)),  yrel_cache{i}(righttroughs_ind(i)), 'm') % trough to right of max peak in magenta
    hold off
    
    % Plot vertical lines denoting max average ABR "chunk" bounds at
    % troughs
    if i==A_length
        hold on
        line([x(chunk1_start) x(chunk1_start)], ylim, 'LineWidth', 1, 'Color', 'k'); % vertical line for chunk1 start
        line([x(chunk2_start) x(chunk2_start)], ylim, 'LineWidth', 1, 'Color', 'k'); % vertical line for chunk2 start
        line([x(chunk3_start) x(chunk3_start)], ylim, 'LineWidth', 1, 'Color', 'k'); % vertical line for chunk2 start
        line([x(chunk3_end) x(chunk3_end)], ylim, 'LineWidth', 1, 'Color', 'k'); % vertical line for chunk2 start
        hold off
    end
    
    % Plot histogram of single-trace inner products
    % Chunk 1
    count2 = count2 + 1;
    axesHandle1(i) = subplot(A_length, 4, count2); 
    histogram(dist_innerprod1{i}, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'edgecolor', 'none', 'FaceColor', 'c')
    xlabel('Inner product (a.u.)')
    ylabel('Prob')
    title(['Input A=', num2str(A_csv(i)), ' dB SPL'])
    
    % plot mean as blue line
    hold on
    signal_mean = mean(dist_innerprod1{i});
    line([signal_mean signal_mean], ylim, 'LineWidth', 1, 'Color', 'b'); % vertical line for cutoff
    hold off
    
    % Chunk 2
    count2 = count2 + 1;
    axesHandle2(i) = subplot(A_length, 4, count2); 
    histogram(dist_innerprod2{i}, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'edgecolor', 'none', 'FaceColor', 'g')
    xlabel('Inner product (a.u.)')
    ylabel('Prob')
    title(['Input A=', num2str(A_csv(i)), ' dB SPL'])
    
    % plot mean as blue line
    hold on
    signal_mean = mean(dist_innerprod2{i});
    line([signal_mean signal_mean], ylim, 'LineWidth', 1, 'Color', 'b'); % vertical line for cutoff
    hold off
    
    % Chunk 3
    count2 = count2 + 1;
    axesHandle3(i) = subplot(A_length, 4, count2); 
    histogram(dist_innerprod3{i}, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'edgecolor', 'none', 'FaceColor', 'm')
    xlabel('Inner product (a.u.)')
    ylabel('Prob')
    title(['Input A=', num2str(A_csv(i)), ' dB SPL'])
    
    % plot mean as blue line
    hold on
    signal_mean = mean(dist_innerprod3{i});
    line([signal_mean signal_mean], ylim, 'LineWidth', 1, 'Color', 'b'); % vertical line for cutoff
    hold off
end
% Same y axes for average ABR plots
same_yaxes(AxesHandles_4)
% Same x and y axes for all histograms
same_yaxes(axesHandle1)
same_xaxes(axesHandle1)
% Same x and y axes for all histograms
same_yaxes(axesHandle2)
same_xaxes(axesHandle2)
% Same x and y axes for all histograms
same_yaxes(axesHandle3)
same_xaxes(axesHandle3)

%% 5-23-19: Calculate p-value for each dB SPL level, using Wilcoxin sign-rank test or K-S test
% compared with 0 dB SPL inner product distribution
% Compute for all chunks!
alldist_innerprod = cell(3,1);
alldist_innerprod{1} = dist_innerprod1;
alldist_innerprod{2} = dist_innerprod2;
alldist_innerprod{3} = dist_innerprod3;
chunk_names = {'Stimulus-1st trough', '1st-2nd trough', '2nd-3rd trough'};

for j = 1:length(alldist_innerprod)
    [p_val, signrank_stat, p_val_KS, ks_score] = innerprod2pval(alldist_innerprod{j}); % 5-26-2019 implementation

    figure('DefaultAxesFontSize', 20)
    % Plot p-value vs dB SPL level: Wilcoxin sign rank test
    subplot(2, 2, 1); 
    plot(A_csv, p_val, '-o')
    title([chunk_names{j}, ': Wilcoxin sign rank for inner products vs 0 dB SPL'])
    xlabel('Amplitude (dB SPL)')
    ylabel('p-value')  
    % CUSTOM add coordinate of threshold point
    THRESH_INDEX = 5; % index of A_csv value that yields p_val < 0.05; from looking at plot
    x_text = A_csv(THRESH_INDEX);
    y_text = p_val(THRESH_INDEX);
    text(x_text - 8, y_text, ['(' num2str(x_text) ', ' num2str(y_text, 2) ')'])

    % Plot Wilcoxign sign rank z-stat vs dB SPL level: Wilcoxin sign rank test
    subplot(2, 2, 2); 
    plot(A_csv, signrank_stat, '-o')
    title([chunk_names{j}, ': Wilcoxin sign rank for inner products vs 0 dB SPL'])
    xlabel('Amplitude (dB SPL)')
    ylabel('Sign rank test statistic')  
    % CUSTOM add coordinate of threshold point
    THRESH_INDEX = 5; % index of A_csv value that yields p_val < 0.05; from looking at plot
    x_text = A_csv(THRESH_INDEX);
    y_text = signrank_stat(THRESH_INDEX);
    text(x_text - 8, y_text, ['(' num2str(x_text) ', ' num2str(y_text, 2) ')'])

    % Plot p-value vs dB SPL level: K-S test
    subplot(2, 2, 3); 
    plot(A_csv, p_val_KS, '-o')
    title([chunk_names{j}, ': Kolmogorov-Smirnov test for inner products vs 0 dB SPL'])
    xlabel('Amplitude (dB SPL)')
    ylabel('p-value')  
    % CUSTOM add coordinate of threshold point
    THRESH_INDEX = 5; % index of A_csv value that yields p_val < 0.05; from looking at plot
    x_text = A_csv(THRESH_INDEX);
    y_text = p_val_KS(THRESH_INDEX);
    text(x_text - 8, y_text, ['(' num2str(x_text) ', ' num2str(y_text, 2) ')'])

    % Plot KS_score vs dB SPL level: K-S test
    subplot(2, 2, 4); 
    plot(A_csv, ks_score, '-o')
    title([chunk_names{j}, ': Kolmogorov-Smirnov test for inner products vs 0 dB SPL'])
    xlabel('Amplitude (dB SPL)')
    ylabel('K-S score')  
    % CUSTOM add coordinate of threshold point
    THRESH_INDEX = 5; % index of A_csv value that yields p_val < 0.05; from looking at plot
    x_text = A_csv(THRESH_INDEX);
    y_text = ks_score(THRESH_INDEX);
    text(x_text - 8, y_text, ['(' num2str(x_text) ', ' num2str(y_text, 2) ')'])
end


%% 5-31-19 Calculate xcov lags using chunk around maximum peak
% avg_lag_xcovExtrap2 = lags_xcov(X_csv, dt, 0.5, A_csv, true, chunk2);
avg_lag_xcovExtrap2 = lags_xcov_localmax(X_csv, dt, 0.5, A_csv, true, chunk2);
dist_innerprod = analyze_innerprod_ABR(X_csv, round(avg_lag_xcovExtrap2/dt));
[p_val, ranksum_stat, p_val_KS, ks_stat] = innerprod2pval(dist_innerprod);
visualize_innerprod_stats(p_val, ranksum_stat, p_val_KS, ks_stat, A_csv)

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
    
%     % 5-25-19: PEAK LAGS: Plot peaks and troughs as overlay on relative average ABR
%     hold on
%     scatter(x(peak_cache{i}), yrel_cache{i}(peak_cache{i}), 'b') % blue peaks
%     scatter(x(trough_cache{i}), yrel_cache{i}(trough_cache{i}), 'r') % red troughs
%     scatter(x(maxpeaks_ind(i)),  yrel_cache{i}(maxpeaks_ind(i)), 'g')
%     hold off

    % 5-25-19: XCOV LAGS: plot max peak lag as black line
    hold on
    peaklag_i = x(maxpeaks_ind(A_length))+avg_lag_xcovExtrap2(i);
    line([peaklag_i peaklag_i], ylim, 'LineWidth', 1, 'Color', 'k'); % vertical line for cutoff
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