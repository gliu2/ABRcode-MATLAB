%% jasaFigures_6_8_19.m
% 
% Create figures for JASA paper on inner product calculation of ABR
% single-trace ensemble threshold
%
% Run this script after "import_ABRcsv_folder.m", which loads Noor's
% single-trace ABR data.
%
% Dependencies: analyze_v2_ABR.m, same_yaxes.m, PTDetect.m, vp2p_abr.m,
%               
% Last edit: 5/26/2019
%
% Author: George Liu

SAMPLES = size(X_csv{1}, 1);
SAMPLING_RATE = 200000; % Hz
dt = 1/SAMPLING_RATE * 1000; % ms
T_total = (SAMPLES-1)*dt; % ms
x = 0:dt:T_total;

%% FIGURE 1
%% Plot averaged ABR waveforms -- the current standard for visually determining ABR threshold

RMS2_averaged_trace = zeros(A_length, 1);
figure('DefaultAxesFontSize', 20)
AxesHandles_1 = zeros(A_length, 1);
for i = 1:A_length
    y = mean(X_csv{i}, 2);
    
    AxesHandles_1(i) = subplot(ceil(A_length/2),2,i);
    plot(x, y)
    title(['Input A=', num2str(A_csv(i)), ' dB SPL'])
    xlabel('Time (ms)')
    ylabel('Voltage (V)')
    
    % Append RMS squared
    RMS2_averaged_trace(i, 1) = analyze_v2_ABR(y);
end
same_yaxes(AxesHandles_1)
    
% subplot(2,2,2)       
% plot(x,y2)        
% title('Subplot 2')
% subplot(2,2,3)      
% plot(x,y3)          
% title('Subplot 3')
% subplot(2,2,4)       
% plot(x,y4)
% title('Subplot 4')

%% 5-13-19: Plot averaged ABR waveforms on a relative scale -- "e.g. individual wave was normalized to the maximal wave within each ABR
% pattern", per Zhou et al. 2006 "Auditory brainstem responses in 10 inbred
% strains of mice" 

figure('DefaultAxesFontSize', 20)
AxesHandles_1 = zeros(A_length, 1);
for i = 1:A_length
    y = mean(X_csv{i}, 2);
    y_rel = y/max(abs(y)); % y vector is scaled to max 1
    
    AxesHandles_1(i) = subplot(ceil(A_length/2),2,i);
    plot(x, y_rel)
    title(['Input A=', num2str(A_csv(i)), ' dB SPL'])
    xlabel('Time (ms)')
    ylabel('Relative voltage (a.u.)')
end
same_yaxes(AxesHandles_1)
    

%% 5-25-19: Threshold method (Oghalai/Ricci labs): Find dB level when max peak-peak amplitude 
% in average ABR > 4 (or 5*) standard deviations above noise floor. 
% References: 
% (1) Wenzel et al. "Laser-induced collagen remodeling and deposition 
% within the basilar membrane of the mouse cochlea". J Biomed Opt (2007).
% (2) Xia et al. "Deficient forward transduction and enhanced reverse transduction
% in the alpha tectorin C1509G human hearing loss mutation", Disease Models
% and Mechanisms (2010). 
% (3)* Cho et al. "Mechanisms of Hearing Loss after Blast Injury to the
% Ear". PLOS ONE (2013).
% (4)* Becker et al. "The presynaptic ribbon maintains vesicle populations 
% at the hair cell afferent fiber synapse". eLIFE (2018).

% Identify peaks and troughs
PEAK_THRESH = 0.1;
yave_cache = cell(A_length, 1);
yrel_cache = cell(A_length, 1);
peak_cache = cell(A_length, 1);
trough_cache = cell(A_length, 1);

figure('DefaultAxesFontSize', 20)
AxesHandles_1 = zeros(A_length, 1);
for i = 1:A_length
    y = mean(X_csv{i}, 2); % average ABR; volts (V)
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
    yave_cache{i} = y;
    yrel_cache{i} = y_rel;
    peak_cache{i} = P;
    trough_cache{i} = T;
end
same_yaxes(AxesHandles_1)

%Find maximum peak-to-peak amplitude in average ABR at each dB level
max_p2p = zeros(A_length, 1);
for i=1:A_length
    max_p2p(i) = vp2p_abr(mean(X_csv{i}, 2), PEAK_THRESH);
end


%% 5-25-19: Plot maximum peak-to-peak amplitude vs dB level
figure('DefaultAxesFontSize', 20)
RMS_averaged_trace_uV = RMS_averaged_trace*10^6; % change standard deviation (RMS) units V -> uV y-axis units
y_uV = max_p2p*10^6; % in case want to change V -> uV y-axis units
NUM_STD_ABOVENOISE = 4; % 4 or 5
thresh_up = NUM_STD_ABOVENOISE*RMS_averaged_trace_uV(1);
thresh_1 = thresh_up + y_uV(1);
plot(A_csv, y_uV, '-o');
my_yline = yline(thresh_1, '--', ['THRESHOLD = ', num2str(NUM_STD_ABOVENOISE), '(RMS)_{0}+V_0^{p-p} = (', num2str(thresh_up, 3), ' + ', num2str(y_uV(1), 3), ')\muV = ', num2str(thresh_1), ' \muV']);
title(['Peak-to-peak amplitude of average ABR'], 'FontSize', 32)
xlabel('Stimulus level (dB SPL)', 'FontSize', 32)
ylabel('Max peak-to-peak amplitude (\muV)', 'FontSize', 32)
set(my_yline, 'FontSize', 18)
% set(gca,'FontSize',20)

% Identify threshold and mark vertical line
isrms_above = y_uV > thresh_1;
if any(y_uV)
    Athresh_ind = find(isrms_above, 1) - 1;
    A_thresh = A_csv(Athresh_ind);
    
    % Linear interpolate to identify exact amplitude corresponding to
    % 3*RMS_0
    dify = y_uV(Athresh_ind + 1) - y_uV(Athresh_ind);
    difx = A_csv(Athresh_ind + 1) - A_thresh; % should be 5 dB or whatever amplitude spacing is
    mslope = dify/difx;
    deltay = thresh_1 - y_uV(Athresh_ind);
    A_tresh_exact = A_thresh + deltay/mslope;
    
    % Draw vertical line at exact threshold amplitude 
    hold on
    line([A_tresh_exact A_tresh_exact], [0, thresh_1], 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--'); % vertical line at exact threshold
    set(gca, 'XTick', sort([round(A_tresh_exact, 0), get(gca, 'XTick')]));
    hold off
else
    disp('Warning: No ABR RMS threshold!')
end


%% 6-3-19: Plot distribution of single-trace Vp-p and innerprods for each dB level
PEAK_THRESH = 0.10;

max_p2p_alldb = cell(A_length, 1);
for i = 1:A_length
%     max_p2p_alldb{i} = vp2p_abr(X_csv{i}, PEAK_THRESH);
    max_p2p_alldb{i} = vp2p_abr_sp(X_csv{i}, PEAK_THRESH); 
end

avg_lag_xcovExtrap = lags_xcov(X_csv, dt, 0.5, A_csv);
dist_innerprod = analyze_innerprod_ABR(X_csv, round(avg_lag_xcovExtrap/dt), 'coeff');

% Plot averaged ABR trace and distribution of single-trace signal
% components per dB SPL level
count2 = 0;
figure('DefaultAxesFontSize', 16)
AxesHandles_3 = zeros(A_length, 1);
AxesHandles_4 = zeros(A_length, 1);
axesHandle = zeros(A_length, 1);
axesHandle2 = zeros(A_length, 1);
for i = 1:A_length
    % Plot averaged ABR trace in first column - absolute scale
    count2 = count2 + 1;
    AxesHandles_3(i) = subplot(A_length, 4, count2);  
    y = mean(X_csv{i}, 2);
    plot(x, y*10^6)
    title(['Average ABR, Input A=', num2str(A_csv(i)), ' dB SPL'])
    xlabel('Time (ms)')
    ylabel('\muV')
    
    % Plot averaged ABR trace in second column - relative scale
    count2 = count2 + 1;
    AxesHandles_4(i) = subplot(A_length, 4, count2);  
    plot(x, y/max(abs(y)))
    title(['Average ABR, Input A=', num2str(A_csv(i)), ' dB SPL'])
    xlabel('Time (ms)')
    ylabel('Rel volt')
    
%     % 5-25-19: XCOV LAGS: plot max peak lag as black line
%     hold on
%     peaklag_i = x(maxpeaks_ind(A_length))+this_lag(i);
%     line([peaklag_i peaklag_i], ylim, 'LineWidth', 1, 'Color', 'k'); % vertical line for cutoff
%     hold off
    
    % Plot histogram of single-trace peak-peak voltages
    count2 = count2 + 1;
    axesHandle(i) = subplot(A_length, 4, count2); 
    histogram(max_p2p_alldb{i}, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'edgecolor', 'r', 'DisplayStyle','stairs')
    hold on
    histogram(max_p2p_alldb{1}, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'edgecolor', 'b', 'DisplayStyle','stairs')
    hold off
    xlabel('Time-aligned peak-peak voltage (V)')
    ylabel('Prob')
    title(['Input A=', num2str(A_csv(i)), ' dB SPL'])
%     xlim([0, 1]*10^-5)
    xlim([-0.6, 0.9]*10^-5)

    % plot mean as blue line
    hold on
    signal_mean = mean(max_p2p_alldb{i});
    line([signal_mean signal_mean], ylim, 'LineWidth', 1, 'Color', 'r'); % vertical line for cutoff
    signal_mean0 = mean(max_p2p_alldb{1});
    line([signal_mean0 signal_mean0], ylim, 'LineWidth', 1, 'Color', 'b'); % vertical line for cutoff
    hold off

    % Plot histogram of single-trace inner products
    count2 = count2 + 1;
    axesHandle2(i) = subplot(A_length, 4, count2); 
    histogram(dist_innerprod{i}, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'edgecolor', 'r', 'DisplayStyle','stairs')
    hold on
    histogram(dist_innerprod{1}, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'edgecolor', 'b', 'DisplayStyle','stairs')
    hold off
    xlabel('Inner product (a.u.)')
    ylabel('Prob')
    title(['Input A=', num2str(A_csv(i)), ' dB SPL'])

    % plot mean as blue line
    hold on
    signal_mean = mean(dist_innerprod{i});
    line([signal_mean signal_mean], ylim, 'LineWidth', 1, 'Color', 'r'); % vertical line for cutoff
    signal_mean0 = mean(dist_innerprod{1});
    line([signal_mean0 signal_mean0], ylim, 'LineWidth', 1, 'Color', 'b'); % vertical line for cutoff
    hold off
end
% Same y axes for average ABR plots
same_yaxes(AxesHandles_3)
same_yaxes(AxesHandles_4)
% Same x and y axes for all histograms
same_yaxes(axesHandle)
same_xaxes(axesHandle)

same_yaxes(axesHandle2)
same_xaxes(axesHandle2)

%% 6-4-19: Psychometric curves for Vp-p and innerproducts of single trace distributions
FP_RATE = 0.05;

dist_allfeatures = {max_p2p_alldb, dist_innerprod};
names_allfeatures = {'V_{p-p}^{mp}', 'innerprod'};

% Vp-p
cut_fp05 = zeros(length(dist_allfeatures), 1);
for j=1:length(dist_allfeatures)
    dist_cell = dist_allfeatures{j};
    noise_dist = dist_cell{1};
    cut_fp05(j) = prctile(noise_dist, (1-FP_RATE)*100); % 5% false positive cutoff

    % Hit rates vs dB
    hits = zeros(A_length, 1);
    hit_rate = zeros(A_length, 1);
    for i = 1:A_length
        this_dist = dist_cell{i};
        hits(i) = sum(this_dist > cut_fp05(j));
%         hits(i) = sum(abs(this_dist) > cut_fp05(j)); % 2-sided tail
        hit_rate(i) = hits(i) / numel(this_dist);
    end

    % Plot 
    figure('DefaultAxesFontSize', 20)
    plot(A_csv, hit_rate, '-o')
    title(['Psychometric curve for ', names_allfeatures{j}, ', FPR=', num2str(FP_RATE)], 'FontSize', 32)
    xlabel('Stimulus level (dB SPL)', 'FontSize', 32)
    ylabel('Hit rate', 'FontSize', 32)
    xline(A_tresh_exact);
    set(gca, 'XTick', sort([round(A_tresh_exact, 0), get(gca, 'XTick')]));

    THRESH_PSYCH = 0.50;
    isrms_above = hit_rate > THRESH_PSYCH;
    if any(hit_rate)
        Athresh_ind3 = find(isrms_above, 1) - 1;
        A_thresh3 = A_csv(Athresh_ind3);

        % Linear interpolate to identify exact amplitude corresponding to
        % 3*RMS_0
        dify = hit_rate(Athresh_ind3 + 1) - hit_rate(Athresh_ind3);
        difx = A_csv(Athresh_ind3 + 1) - A_thresh3; % should be 5 dB or whatever amplitude spacing is
        mslope = dify/difx;
        deltay = THRESH_PSYCH - hit_rate(Athresh_ind3);
        A_thresh_exact3 = A_thresh3 + deltay/mslope;

        % Draw horizontal line at FPR threshold, vertical line at exact threshold amplitude 
        hold on
        my_yline = yline(THRESH_PSYCH, '--', ['THRESHOLD = ', num2str(THRESH_PSYCH, 3)]);
        set(my_yline, 'FontSize', 18)
        line([A_thresh_exact3 A_thresh_exact3], [0, THRESH_PSYCH], 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--'); % vertical line at exact threshold
        set(gca, 'XTick', sort(unique([round(A_thresh_exact3, 0), get(gca, 'XTick')])));
        hold off
    else
        disp('Warning: No ABR RMS threshold!')
    end
end

%% 6-5-19: Plot Wilcoxin rank-sum test p-values for each metric
for j=1:length(dist_allfeatures)
    dist_cell = dist_allfeatures{j};
    
    % Get Wilcoxin sum-rank p-values
    [p_val, ranksum_stat, p_val_KS, ks_stat] = innerprod2pval(dist_cell);
    
    % Plot 
    figure('DefaultAxesFontSize', 20)
    plot(A_csv, p_val, '-o')
    title(['P-value curve for ', names_allfeatures{j}, ', rank-sum'], 'FontSize', 32)
    xlabel('Stimulus level (dB SPL)', 'FontSize', 32)
    ylabel('P-value', 'FontSize', 32)
    xline(A_tresh_exact);
    set(gca, 'XTick', sort([round(A_tresh_exact, 0), get(gca, 'XTick')]));
    
    % Plot 
    figure('DefaultAxesFontSize', 20)
    plot(A_csv, p_val_KS, '-o')
    title(['P-value curve for ', names_allfeatures{j}, ', K-S'], 'FontSize', 32)
    xlabel('Stimulus level (dB SPL)', 'FontSize', 32)
    ylabel('P-value', 'FontSize', 32)
    xline(A_tresh_exact);
    set(gca, 'XTick', sort([round(A_tresh_exact, 0), get(gca, 'XTick')]));
    
end

%% 6-6-19: On peakpeak average ABR vs db, Identify threshold and mark vertical line
thresh_2 = cut_fp05(1)*10^6;

figure('DefaultAxesFontSize', 20)
RMS_averaged_trace_uV = RMS_averaged_trace*10^6; % change standard deviation (RMS) units V -> uV y-axis units
max_p2p_uV = max_p2p*10^6; % in case want to change V -> uV y-axis units
h8 = plot(A_csv, max_p2p_uV, '-o');
title(['Peak-to-peak amplitude of average ABR'], 'FontSize', 32)
xlabel('Stimulus level (dB SPL)', 'FontSize', 32)
ylabel('Max peak-to-peak amplitude (\muV)', 'FontSize', 32)


isrms_above = max_p2p_uV > thresh_2;
if any(isrms_above)
    A_tresh_exact2 = xgiveny(thresh_2, max_p2p_uV, A_csv);
%     Athresh_ind2 = find(isrms_above, 1) - 1;
%     A_thresh2 = A_csv(Athresh_ind2);
%     
%     % Linear interpolate to identify exact amplitude corresponding to
%     % 3*RMS_0
%     dify = max_p2p_uV(Athresh_ind2 + 1) - max_p2p_uV(Athresh_ind2);
%     difx = A_csv(Athresh_ind2 + 1) - A_thresh2; % should be 5 dB or whatever amplitude spacing is
%     mslope = dify/difx;
%     deltay = thresh_2 - max_p2p_uV(Athresh_ind2);
%     A_tresh_exact2 = A_thresh2 + deltay/mslope;
    
    % Draw horizontal line at FPR threshold, vertical line at exact threshold amplitude 
    hold on
    my_yline = yline(thresh_2, '--', ['CUTOFF = ', num2str(thresh_2, 3), ' \muV for FPR ', num2str(FP_RATE)]);
    set(my_yline, 'FontSize', 18)
    line([A_tresh_exact2 A_tresh_exact2], [0, thresh_2], 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--'); % vertical line at exact threshold
    set(gca, 'XTick', sort([round(A_tresh_exact2, 0), get(gca, 'XTick')]));
    hold off
else
    disp('Warning: No ABR RMS threshold!')
end

%% 6-6-19: inner products on average ABR vs dB level
thresh_4 = cut_fp05(2);

X_csv_aveonly = cell(A_length, 1);
for i = 1:A_length
    X_csv_aveonly{i} = mean(X_csv{i}, 2);
end

avg_lag_xcovExtrap_ave = lags_xcov(X_csv_aveonly, dt, 0.5, A_csv);
dist_innerprod_ave = analyze_innerprod_ABR(X_csv_aveonly, round(avg_lag_xcovExtrap_ave/dt), 'coeff');
dist_innerprod_ave = cell2mat(dist_innerprod_ave);
% disp(size(dist_innerprod_ave))

figure('DefaultAxesFontSize', 20)
plot(A_csv, dist_innerprod_ave, '-o');
title(['Inner product of average ABR'], 'FontSize', 32)
xlabel('Stimulus level (dB SPL)', 'FontSize', 32)
ylabel('Inner product (a.u.)', 'FontSize', 32)


isrms_above = dist_innerprod_ave > thresh_4;
if any(isrms_above)
    Athresh_ind4 = find(isrms_above, 1) - 1;
    A_thresh4 = A_csv(Athresh_ind4);
    
    % Linear interpolate to identify exact amplitude corresponding to
    % 3*RMS_0
    dify = dist_innerprod_ave(Athresh_ind4 + 1) - dist_innerprod_ave(Athresh_ind4);
    difx = A_csv(Athresh_ind4 + 1) - A_thresh4; % should be 5 dB or whatever amplitude spacing is
    mslope = dify/difx;
    deltay = thresh_4 - dist_innerprod_ave(Athresh_ind4);
    A_tresh_exact4 = A_thresh4 + deltay/mslope;
    
    % Draw horizontal line at FPR threshold, vertical line at exact threshold amplitude 
    hold on
    my_yline = yline(thresh_4, '--', ['CUTOFF = ', num2str(thresh_4, 3), ' \muV for FPR ', num2str(FP_RATE)]);
    set(my_yline, 'FontSize', 18)
    line([A_tresh_exact4 A_tresh_exact4], [0, thresh_4], 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--'); % vertical line at exact threshold
    set(gca, 'XTick', sort(unique([round(A_tresh_exact4, 0), get(gca, 'XTick')])));
    hold off
else
    disp('Warning: No ABR RMS threshold!')
end

%% 6-7-19: ROC curves for single-trace features
for i = 1:length(dist_allfeatures)
    this_feature = dist_allfeatures{i};
    noise_dist = this_feature{1};
    feature_name = names_allfeatures{i};
    
    
    % Get ROC Curve for feature i, one for each dB level
    my_colors = jet(A_length);
    feature_AUC = zeros(A_length, 1);
    figure('DefaultAxesFontSize', 16)
    hold on
    for j=1:A_length
        A_dist = this_feature{j};
        scores = [noise_dist(:); A_dist(:)];
        labels = [zeros(numel(noise_dist), 1); ones(numel(A_dist), 1)];
        posclass = 1;
        
        [X_roc, Y_roc,T_roc,feature_AUC(j)] = perfcurve(labels,scores,posclass);
                 
        plot(X_roc, Y_roc, 'Color', my_colors(j, :))
       
    end
    legend([num2str(A_csv), repmat(' dB, AUC=', A_length, 1), num2str(feature_AUC)])
    xlabel('False positive rate'); ylabel('True positive rate');
    title(['ROC Curves for ', feature_name])
    hold off
end

%% 6-7-19: FPR vs thresh curves for average trace features

% Cell array of each feature of average ABRs vs dB
dist_avefeatures = {max_p2p_uV*10^-6, dist_innerprod_ave};

for i = 1:length(dist_allfeatures)
    this_feature_ave = dist_avefeatures{i};
    this_feature = dist_allfeatures{i};
    noise_dist = this_feature{1};
    feature_name = names_allfeatures{i};
    
    % Get FPR vs cutoff
    [f_fpr, x_cuts] = ecdf(noise_dist);
    f_fpr = 1 - f_fpr; 
    
    % Get cutoff vs threshold for average trace data - threshold is dB where
    % feature of mean trace > cutoff
    thresh_curve = zeros(size(x_cuts));
    for j = 1:length(x_cuts)
        this_cut = x_cuts(j);
        thresh_curve(j) = xgiveny(this_cut, this_feature_ave, A_csv);
    end
    
    % Plot FPR (y) vs threshold (x)
    figure('DefaultAxesFontSize', 20)
    plot(thresh_curve, f_fpr)
    xlabel('Threshold (dB)'); 
    ylabel('False positive rate');
    title(['FPR vs threshold for ', feature_name])
    
    % Plot FPR (y) vs cutoff (x)
    figure('DefaultAxesFontSize', 20)
    plot(x_cuts, f_fpr)
    xlabel('Cutoff'); 
    ylabel('False positive rate');
    title(['FPR vs cutoff for ', feature_name])
    
    % Plot cutoff (y) vs threshold (x)
    figure('DefaultAxesFontSize', 20)
    plot(thresh_curve, x_cuts)
    xlabel('Threshold (dB)'); 
    ylabel('Cutoff');
    title(['Cutoff vs threshold for ', feature_name])
end