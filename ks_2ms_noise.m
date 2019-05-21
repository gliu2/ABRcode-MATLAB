%% ks_2ms_noise.m
% 
% Analyze Noor's ABR traces with 2 ms noise at beginning.
%
% Run this script after running 'import_ABRcsv_folder.m'
%
% Dependencies: same_yaxes.m
% Last edit: 5/17/2019 - Analyze noise distribution 
%
% Author: George Liu

X95 = X_csv{end};
t_endnoise = find(x==2);
% t_endnoise = size(X95, 1);
noise_only = X95(1:t_endnoise, :);
mm = size(noise_only, 2);

KS_plot = zeros(mm-1, 1);
p_plot = zeros(mm-1, 1);
for ii = 2:mm
    x1 = noise_only(:, ii-1);
    x2 = noise_only(:, ii);
    [h,p,ks2stat] = kstest2(x1, x2);
    
    p_plot(ii-1, 1) = p;
    KS_plot(ii-1, 1) = ks2stat;
end

figure
plot(KS_plot, 'LineWidth', 1)
xlabel('n-1', 'FontSize', 32)
ylabel('KS(Pn, Pn-1)', 'FontSize', 32)
title('Noise KS test, 95 dB SPL', 'FontSize', 32)
% title('Noise KS test, 0 dB SPL', 'FontSize', 32)

figure
plot(p_plot, 'LineWidth', 1)
xlabel('n-1', 'FontSize', 32)
ylabel('P-value KS(Pn, Pn-1)', 'FontSize', 32)
title('Noise KS p-value, 95 dB SPL', 'FontSize', 32)
% title('Noise KS p-value, 0 dB SPL', 'FontSize', 32)

%% 5-16-19: Calculate K-S statistic distribution histogram, for aggregated across all dB levels (2 ms noise inter-repetition comparisons)
allKS_stat = [];
all_p = [];
% Iterate over stimulus levels (dB SPL)
for i = 1:A_length
    this_X = X_csv{i};
    this_noise = this_X(1:t_endnoise, :); % 2ms noise at beginning
%     this_m = size(this_X, 2);
    
%     KS_thisdb = zeros(mm-1, 1);
%     p_plot = zeros(mm-1, 1);
    % Iterate over consecutive single trace pairs
    for j = 2:mm
    x1 = this_noise(:, ii-1);
    x2 = this_noise(:, ii);
    [~, this_p, this_ks2stat] = kstest2(x1, x2);
    allKS_stat = [allKS_stat; this_ks2stat];
    all_p = [all_p; this_p];
    end
end

% Calculate KS stat corresponding to critical p value
P_CRIT = 0.01;
[~, ind] = min(abs(all_p - P_CRIT));
ks_crit = allKS_stat(ind);
% num_below = 

% Plot histogram of KS test statistics
figure('DefaultAxesFontSize', 20)
histogram(allKS_stat, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'FaceColor', 'r')
ylabel('Probability')
xlabel('KS2 test stat')
title(['Histogram of Kolmogorov-Smirnov Statistic at all dB, p_{crit}=', num2str(P_CRIT), ' Kstat_{crit}=', num2str(ks_crit)])

% critical p-value / KS stat as blue line
hold on
line([ks_crit ks_crit], ylim, 'LineWidth', 1, 'Color', 'b'); % vertical line for p_crit
% set(gca, 'XTick', sort([round(ks_crit, 1), get(gca, 'XTick')]));
hold off

%% 5-17-19: Measure correlation time
this_X = X_csv{1};
cc_time = zeros(t_endnoise, 1);
TIME_PT = 200;
for i = 1:t_endnoise
    cc = corrcoef(this_X(TIME_PT,:), this_X(i, :));
    cc_time(i) = cc(1,2);
end

figure('DefaultAxesFontSize', 20)
plot(x(1:t_endnoise)-TIME_PT*dt, cc_time)
xlabel('\Deltat (ms)', 'FontSize', 32)
ylabel('Correlation coefficient', 'FontSize', 32)
title(['Corrcoef of noise for single points \Deltat apart, ', num2str(A_csv(1)), ' dB SPL'], 'FontSize', 32)

%% 5-19-19: Measure autocorrelation of noise traces at 0 dB level
this_m = m(1); % 0 dB single traces
% this_X = X_csv{1}; % SAMPLES x m_traces matrix
this_X = X_csv{1}(1:t_endnoise, :); % NOISE_SAMPLES x m_traces matrix
this_acf = zeros(2*size(this_X,1)-1, this_m);
for i = 1:this_m
    trace_i = this_X(:, i);
    [this_acf(:,i), lags] = xcorr(trace_i); % autocorrelation
end

% Plot random single-trace autocorrelation functions
NUM_SUBPLOT = 6;
rnd_traces = sort(randperm(this_m, NUM_SUBPLOT^2));
% rnd_traces = ceil(rand(NUM_SUBPLOT)*this_m);
% rnd_traces = reshape(sort(rnd_traces(:)), NUM_SUBPLOT, NUM_SUBPLOT)'; % sort in increasing order, row by row
figure
AxesHandles_1 = zeros(NUM_SUBPLOT^2, 1);
count = 0;
for i = 1:NUM_SUBPLOT
    for j = 1:NUM_SUBPLOT
        count = count + 1;
        this_trace = rnd_traces(count);
        AxesHandles_1(count) = subplot(NUM_SUBPLOT, NUM_SUBPLOT, count);
        plot(lags*dt, this_acf(:, this_trace))
        xlabel('Lag (ms)', 'FontSize', 12)
        ylabel('ACF (V^2)', 'FontSize', 12)
        title(['Trace ', num2str(this_trace), ', ', num2str(A_csv(1)), ' dB'])
    end
end
same_yaxes(AxesHandles_1)

% Plot average ACF
figure('DefaultAxesFontSize', 20)
plot(lags*dt, mean(this_acf,2))
xlabel('Lag (ms)', 'FontSize', 32)
ylabel('ACF ave (V^2)', 'FontSize', 32)
title(['Autocorrelation average, ', num2str(A_csv(1)), ' dB SPL'], 'FontSize', 32)

% Match y-axes of single and average ACF plots
allYLim = get(AxesHandles_1, {'YLim'});
allYLim = cat(2, allYLim{:});
ylim([min(allYLim), max(allYLim)]);

%% 5-19-19: PSD (periodogram) of single and average ABR traces at 0 dB
this_m = m(1); % 0 dB single traces
this_X = X_csv{1}; % SAMPLES x m_traces matrix
% this_X = X_csv{1}(1:t_endnoise, :); % NOISE_SAMPLES x m_traces matrix
nfft = length(x);
[pxx, f] = periodogram(this_X, [], nfft, SAMPLING_RATE); % pxx - [] x m_traces matrix, f -  [0,fs/2] (Hz)

% Plot random single-trace autocorrelation functions
NUM_SUBPLOT = 6;
rnd_traces = sort(randperm(this_m, NUM_SUBPLOT^2));
figure
AxesHandles_1 = zeros(NUM_SUBPLOT^2, 1);
count = 0;
for i = 1:NUM_SUBPLOT
    for j = 1:NUM_SUBPLOT
        count = count + 1;
        this_trace = rnd_traces(count);
        AxesHandles_1(count) = subplot(NUM_SUBPLOT, NUM_SUBPLOT, count);
        plot(f, pxx(:, this_trace))
        xlabel('Freq (Hz)', 'FontSize', 12)
        ylabel('PSD (V^2 Hz^{-1})', 'FontSize', 12)
        title(['Trace ', num2str(this_trace), ', ', num2str(A_csv(1)), ' dB'])
        xlim([0, 5000])
    end
end
same_yaxes(AxesHandles_1)

% Plot average ACF
figure('DefaultAxesFontSize', 20)
errorbar(f, mean(pxx, 2), std(pxx, 0, 2)/sqrt(this_m))
% plot(f, mean(pxx,2))
xlabel('Freq (Hz)', 'FontSize', 32)
ylabel('PSD (V^2 Hz^{-1})', 'FontSize', 32)
title(['PSD average, ', num2str(A_csv(1)), ' dB SPL'], 'FontSize', 32)
xlim([0, 5000])

% Match y-axes of single and average ACF plots
allYLim = get(AxesHandles_1, {'YLim'});
allYLim = cat(2, allYLim{:});
ylim([min(allYLim), max(allYLim)]);

%% Subplots for individual trace ABRs
figure % plot 10 x 5 subplots. A_length = 10
num_examples = 5;
for j = 1:num_examples
    y2 = X95(:, j); % loop over first <num_examples> single traces
    subplot(num_examples, 1, j)  
    plot(x, y2)
    xlabel('Time (ms)')
    ylabel('Voltage (V)')
    title(['Input A=', num2str(95), ' dB SPL'])
end

%% Subplots for individual trace ABRs histograms
figure % plot 10 x 5 subplots. A_length = 10
num_examples = 5;
for j = 1:num_examples
    y2 = noise_only(:, j); % loop over first <num_examples> single traces
    subplot(5, 1, j)  
    histogram(y2, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'FaceColor', 'r')
    ylabel('Probability')
    xlabel('Voltage (V)')
    title(['Input A=', num2str(95), ' dB SPL, trace ', num2str(j)])
end

%% Find threshold versus time based on single-point analysis
FP_RATE = 0.20;
t_endnoise = find(x==2);
    
% TIME = 0.04;
% t_TIME = find(x==4);

num_timepts = size(X_csv{end}, 1) - size(noise_only, 1);
prob_hit_all = zeros(num_timepts, 20);
for aa = 1:20
    disp(['Working on ', num2str(aa), ' out of ', num2str(20), '...'])
    this_X = X_csv{aa}; % one level at max 95 dB SPL
    t_first = t_endnoise + 1;
    t_last = size(this_X, 1);
    num_timepts = length(t_first:t_last);
    % prob_hit = zeros(num_timepts, 1);

    % % For each amplitude, calculate hit rate
    % prob_hit = zeros(A_length, 1);
    % for aa = 1:A_length
    %     this_X  = X_csv{aa};

    % For each single trace:
    hit_or_miss = zeros(num_timepts, mm);
    for iii = 1:mm
        % Find thresholds for false positive rate
        noise_2ms = this_X(1:t_endnoise, iii); 
        volt_thresholds = prctile(noise_2ms, [FP_RATE, (1-FP_RATE)]*100);
        th_high = volt_thresholds(2);
        th_low = volt_thresholds(1);

        % For each time point:
        for tt = t_first:t_last

            % Determine if ABR signal is hit or miss
    %         signal_TIME = this_X(t_TIME, iii); % single time point ABR response value
            signal_TIME = this_X(tt, iii); % single time point ABR response value
            hit_or_miss(tt - t_endnoise, iii) = signal_TIME > th_high || signal_TIME < th_low;
        end
    end
    %     prob_hit(aa) = mean(hit_or_miss);
    prob_hit = mean(hit_or_miss, 2);
    prob_hit_all(:, aa) = prob_hit;
    
    % % Find SPL level that gives 50% hits at voltage detection threshold above 
    % A_thresh = ZeroGL(A_csv, prob_hit-0.5);

    % figure
    % plot(A_csv, prob_hit, 'LineWidth', 1)
    % xlabel('Amplitude (dB SPL)', 'FontSize', 32)
    % ylabel('Probability of hit', 'FontSize', 32)
    % title('Psychometric curve @ 4ms (2ms after sound pip)', 'FontSize', 32)

    figure
    time_pts = x(t_first:t_last);
    plot(time_pts, prob_hit, 'LineWidth', 1)
    xlabel('Time (ms)', 'FontSize', 32)
    ylabel('Probability of hit', 'FontSize', 32)
%     title('Hit probability @ 95 dB SPL', 'FontSize', 32)
    title(['Hit probability @ ', num2str(A_csv(aa)), ' dB SPL'], 'FontSize', 32)
end

figure
xxx = repmat(x(t_first:t_last)',1,A_length);
yyy = repmat(A_csv', num_timepts, 1);

h1 = surf(xxx, yyy, prob_hit_all);
set(h1,'LineStyle','none')
xlabel('Time (ms)')
ylabel('dB SPL')
zlabel('Probability hit')
title('Hit probability @ single time point')
colorbar

%% ALTERNATIVE 4-18-19: noise dist over all traces (Find threshold versus time based on single-point analysis)
FP_RATE = 0.05;
t_endnoise = find(x==2);
    
% TIME = 0.04;
% t_TIME = find(x==4);

prob_hit_all = zeros(num_timepts, 20);
for aa = 1:20
    disp(['Working on ', num2str(aa), ' out of ', num2str(20), '...'])
    this_X = X_csv{aa}; % one level at max 95 dB SPL
    
    % 4-22-19 combine opposite polarity pairs
    X_odd = this_X(:, 1:2:end); % odd rows
    X_even = this_X(:, 2:2:end); % even rows
    X_addpolar = (X_odd + X_even)/2;
    this_X = X_addpolar;
    
    t_first = t_endnoise + 1;
    t_last = size(this_X, 1);
    num_timepts = length(t_first:t_last);
    % prob_hit = zeros(num_timepts, 1);

    % % For each amplitude, calculate hit rate
    % prob_hit = zeros(A_length, 1);
    % for aa = 1:A_length
    %     this_X  = X_csv{aa};

    % Find thresholds for false positive rate
    noise_2ms = this_X(1:t_endnoise, :); 
    volt_thresholds = prctile(noise_2ms(:), [FP_RATE, (1-FP_RATE)]*100);
    th_high = volt_thresholds(2);
    th_low = volt_thresholds(1);
    
    % For each single trace:
    mm = size(this_X, 2);
    hit_or_miss = zeros(num_timepts, mm);
    for iii = 1:mm

        % For each time point:
        for tt = t_first:t_last

            % Determine if ABR signal is hit or miss
    %         signal_TIME = this_X(t_TIME, iii); % single time point ABR response value
            signal_TIME = this_X(tt, iii); % single time point ABR response value
%             hit_or_miss(tt - t_endnoise, iii) = signal_TIME > th_high || signal_TIME < th_low; % two-tailed cutoff
            hit_or_miss(tt - t_endnoise, iii) = signal_TIME > th_high; % one-tailed cutoff
        end
    end
    %     prob_hit(aa) = mean(hit_or_miss);
    prob_hit = mean(hit_or_miss, 2);
    prob_hit_all(:, aa) = prob_hit;
    
    % % Find SPL level that gives 50% hits at voltage detection threshold above 
    % A_thresh = ZeroGL(A_csv, prob_hit-0.5);

    % figure
    % plot(A_csv, prob_hit, 'LineWidth', 1)
    % xlabel('Amplitude (dB SPL)', 'FontSize', 32)
    % ylabel('Probability of hit', 'FontSize', 32)
    % title('Psychometric curve @ 4ms (2ms after sound pip)', 'FontSize', 32)

    figure
    time_pts = x(t_first:t_last);
    plot(time_pts, prob_hit, 'LineWidth', 1)
    xlabel('Time (ms)', 'FontSize', 32)
    ylabel('Probability of hit', 'FontSize', 32)
%     title('Hit probability @ 95 dB SPL', 'FontSize', 32)
    title(['Hit probability @ ', num2str(A_csv(aa)), ' dB SPL'], 'FontSize', 32)
end

figure
xxx = repmat(x(t_first:t_last)',1,A_length);
yyy = repmat(A_csv', num_timepts, 1);

h1 = surf(xxx, yyy, prob_hit_all);
set(h1,'LineStyle','none')
xlabel('Time (ms)')
ylabel('dB SPL')
zlabel('Probability hit')
title('Hit probability @ single time point')
colorbar

%% 4-22-19: combine opposite polarity traces
X95 = X_csv{end}; % 1952 x 514 double

% combine consecutive pairs of traces, which TDT ABR system uses tone pips
% of opposite polarities to generate
X95_odd = X95(:, 1:2:end); % odd rows
X95_even = X95(:, 2:2:end); % even rows
X95_addpolar = (X95_odd + X95_even)/2;
X95_minuspolar = (X95_odd - X95_even)/2;

% Plot single traces from original and new ABR responses
figure % plot 10 x 5 subplots. A_length = 10
num_examples = 20;
AxesHandles_1 = zeros(num_examples,1);
for j = 1:num_examples
    y2 = X95(:, j); % loop over first <num_examples> single traces
    AxesHandles_1(j) = subplot(num_examples/2, 2, j);
    plot(x, y2)
    xlabel('Time (ms)')
    ylabel('Voltage (V)')
    title(['Original trace ', num2str(j), ', ', num2str(95), ' dB SPL'])
end
same_yaxes(AxesHandles_1)

figure % plot 10 x 5 subplots. A_length = 10
num_examples = 10;
AxesHandles_2 = zeros(num_examples,1);
for j = 1:num_examples
    y2 = X95_addpolar(:, j); % loop over first <num_examples> single traces
    AxesHandles_2(j) = subplot(num_examples, 1, j);
    plot(x, y2)
    xlabel('Time (ms)')
    ylabel('Voltage (V)')
    title(['Add polarities, ', num2str(95), ' dB SPL'])
end
same_yaxes(AxesHandles_2)

figure % plot 10 x 5 subplots. A_length = 10
num_examples = 10;
AxesHandles_3 = zeros(num_examples,1);
for j = 1:num_examples
    y2 = X95_minuspolar(:, j); % loop over first <num_examples> single traces
    AxesHandles_3(j) = subplot(num_examples, 1, j); 
    plot(x, y2)
    xlabel('Time (ms)')
    ylabel('Voltage (V)')
    title(['Subtract polarities, ', num2str(95), ' dB SPL'])
end
same_yaxes(AxesHandles_3)

%% For entire noise 2 ms distributions, plot variance of noise wrt 
t_endnoise = find(x==2);
    
% TIME = 0.04;
% t_TIME = find(x==4);

for aa = 1:20
    disp(['Working on ', num2str(aa), ' out of ', num2str(20), '...'])
    this_X = X_csv{aa}; % one level at max 95 dB SPL
    noise_2ms = this_X(1:t_endnoise, :); 
    variance_noiseall = var(noise_2ms); % 1x514 vector of variances of 2ms noise per single trace
    
%     % Plot noise variance vs trace number
%     figure
%     plot(variance_noiseall);
%     xlabel('Single trace #', 'FontSize', 28)
%     ylabel('Variance of noise (V^2)', 'FontSize', 28)
%     title(['Noise variance, ', num2str(A_csv(aa)), ' dB SPL'], 'FontSize', 28)

    % Check if variances same for consecutive noise traces, using F-test
    p_plot = zeros(mm-1, 1);
    for ii = 1:mm-1
        x1 = noise_only(:, ii);
        x2 = noise_only(:, ii+1);
        [h,p] = vartest2(x1, x2);

        p_plot(ii, 1) = p;
    end
    figure
    plot(p_plot, 'LineWidth', 1)
    xlabel('n', 'FontSize', 32)
    ylabel('P-value (n, n+1)', 'FontSize', 32)
    title(['Noise F-test p-value, ', num2str(A_csv(aa)), ' dB SPL'], 'FontSize', 32)
end
