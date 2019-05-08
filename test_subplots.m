% run after import_ABRcsv_folder.m
%
% Assesses various "scores" of individual trace ABR for hearing threshold, such as single trace
% amplitude histogram, autocorrelation function, critical time (for 1/e
% decay of exponential fit to autocorrelation function upper envelope), and
% power spectrum density.
%
% We tried to find patterns in single trace data that match the observed
% average trace ABR pattern. But no clear pattern so far.
%
% 3-13-2019
% rewrite for 10 hard-coded csv files
% Dependencies: acf.m, same_yaxes.m
%
% Last edit: 5-7-19 - update make y-axes same
% George S. Liu

SAMPLES = size(X_csv{i}, 1);
SAMPLING_RATE = 200000; % Hz
dt = 1/SAMPLING_RATE * 1000; % ms
T_total = (SAMPLES-1)*dt; % ms
x = 0:dt:T_total;
% y1 = sin(x);
% y2 = cos(x);
% y3= tan(x);
% y4=1./cos(x);

%% Plot averaved ABR waveforms -- the current standard for visually determining ABR threshold

RMS2_averaged_trace = zeros(A_length, 1);
figure
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

%% Subplots for individual trace ABRs
figure % plot 10 x 5 subplots. A_length = 10
num_examples = 5;
% RMS2_signal_traces = [];
% RMS2_noise_traces = [];
RMS2_single_traces = cell(A_length, 1);
AxesHandles_2 = zeros(A_length*num_examples, 1);
for i = 1:A_length
    this_X = X_csv{i}; % 1952 x ? matrix
    for j = 1:num_examples
        y2 = this_X(:, j); % loop over first <num_examples> single traces

        count = (i-1)*num_examples + j;
        AxesHandles_2(count) = subplot(A_length, num_examples, count);
        plot(x, y2)
        title(['Input A=', num2str(A_csv(i)), ' dB SPL'])
        xlabel('Time (ms)')
        ylabel('Voltage (V)')
    end
    
    
    % Concatenate list of single trace RMS^2
    this_RMS2 = analyze_v2_ABR(this_X); % mx1 vector
%     if A_csv(i)==0
%         RMS2_noise_traces = [RMS2_noise_traces; this_RMS2];
%     else
%         RMS2_signal_traces = [RMS2_signal_traces; this_RMS2];
%     end
    RMS2_single_traces{i} = this_RMS2;
end
same_yaxes(AxesHandles_2)

%% Histograms single trace
% create a default color map ranging from red to light pink
length_c = A_length;
red = [1, 0, 0];
blue = [0, 0, 1];
colors_p = [linspace(red(1),blue(1),length_c)', linspace(red(2),blue(2),length_c)', linspace(red(3),blue(3),length_c)'];

% Histograms of Voltage RMS squared
figure
% histogram(RMS2_noise_traces, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'FaceColor', 'r'), xlabel('<V^2>_t (\muV^2)', 'FontSize', 32), ylabel('Prob(<V^2>_t)', 'FontSize', 32), title('ABR histograms of RMS^2', 'FontSize', 40)
hold on
for i = 1:A_length
    this_RMS2 = RMS2_single_traces{i};
    histogram(this_RMS2, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'FaceColor', colors_p(i, :))
end
line([Vd_new Vd_new], [0 0.3], 'LineWidth', 1); % vertical line for cutoff
hold off
box off % remove ticks on top and right borders of plot
% axis tight % makes edges of data flush with left and right borders of plot
legend('Noise', ['Signal (A=', num2str(A_csv(2)'), ')'], ['Signal (A=', num2str(A_csv(3)'), ')'], ['Signal (A=', num2str(A_csv(4)'), ')'], ['Signal (A=', num2str(A_csv(5)'), ')'], ['Signal (A=', num2str(A_csv(6)'), ')'], ['Signal (A=', num2str(A_csv(7)'), ')'], ['Signal (A=', num2str(A_csv(8)'), ')'], ['Signal (A=', num2str(A_csv(9)'), ')'], ['Signal (A=', num2str(A_csv(10)'), ')'], 'location', 'northeast')
legend boxoff % remove box around legend
xlabel('V_{RMS}^2 (V^2)', 'FontSize', 32);
ylabel('Probability', 'FontSize', 32);
title('Single traces', 'FontSize', 32);
set(gca,'FontSize',20)

%% Histograms averaged trace
figure
hold on
yy = ones(A_length, 1);
scatter(RMS2_averaged_trace, yy, 100, colors_p);
% for i = 1:A_length
%     this_RMS2 = RMS2_averaged_trace(i, 1);
%     histogram(this_RMS2, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'FaceColor', colors_p(i, :))
% end
line([Vd_old Vd_old], [0 2], 'LineWidth', 1); % vertical line for cutoff
hold off
box off % remove ticks on top and right borders of plot
% axis tight % makes edges of data flush with left and right borders of plot
legend('RMS^2 (red noise, blue max signal)', ['Vd_{old} (=', num2str(Vd_old), ' nV^2)'], 'location', 'northeast')
legend boxoff % remove box around legend
xlabel('V_{RMS}^2 (nV^2)', 'FontSize', 32);
ylabel('Probability', 'FontSize', 32);
title('Averaged trace', 'FontSize', 32);
ylim([0.9, 1.1])
set(gca,'FontSize',20)

%% Distribution amplitudes at single time point changes with level
DELAY_TONEPIP = 2; % ms when tone pip starts; ADJUST PER EXPERIMENT

i_beforepip = find(x==1); % noisy point before tone pip, REMOVE IF TONE PIP NOT DELAYED 5-8-19
i_noise = find(x==0.5 + DELAY_TONEPIP); % noisy point t < latency peak 1
i_neg_peak1 = find(x==1.745 + DELAY_TONEPIP); % first negative peak/trough
i_inflection_mid = find(x==2.055 + DELAY_TONEPIP); % zero crossing between trough and peak 2/3
i_pos_peak1 = find(x==2.255 + DELAY_TONEPIP); % positive peak wave 2 (3?)

indexes_sp = [i_beforepip, i_noise, i_neg_peak1, i_inflection_mid, i_pos_peak1];
% indexes_sp = [i_noise, i_neg_peak1, i_inflection_mid, i_pos_peak1];
num_sp = numel(indexes_sp);

count2 = 0;
figure
% Loop over signal levels for each time point
AxesHandles_3 = zeros(A_length, 1);
for i = 1:A_length
    this_A = A_csv(i); % signal level
    this_X = X_csv{i}; % 1952 x ? matrix
    
    % Plot averaged ABR trace in first column
    count2 = count2 + 1;
    AxesHandles_3(i) = subplot(A_length, num_sp + 1, count2);  
    y = mean(X_csv{i}, 2);
    plot(x, y)
    title(['Input A=', num2str(A_csv(i)), ' dB SPL'])
    xlabel('Time (ms)')
    ylabel('Voltage (V)')
    hold on
    % Draw red vertical lines at amplitude distributions' time points
    for k = 1:num_sp
        x_sp = x(indexes_sp(k)); % time point for plotting amplitude distribution
%         line([x_sp x_sp], ylim, 'LineWidth', 1, 'Color', 'r'); % vertical line for cutoff

        % If setting ylims same for all subplots, then make vertical red
        % lines all very tall
        this_ylim = ylim;
        line([x_sp x_sp], [-1,1], 'LineWidth', 1, 'Color', 'r'); % vertical line for cutoff
        ylim(this_ylim)
    end
    hold off
    
    % Loop over single time points
    for j = 1:num_sp
        x_sp = x(indexes_sp(j)); % time point for plotting amplitude distribution
        A_distribution = this_X(indexes_sp(j), :);
        
        count2 = count2 + 1;
        subplot(A_length, num_sp + 1, count2)  
        histogram(A_distribution, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'FaceColor', 'r')
        title(['T = ', num2str(x_sp), ' ms'])
        xlabel('Amplitude (V)')
        ylabel('Probability')
        xlim([-3, 3]*10^-6)
        
        % plot mean as blue line
        hold on
        A_mean = mean(A_distribution);
        line([A_mean A_mean], ylim, 'LineWidth', 1, 'Color', 'b'); % vertical line for cutoff
        hold off
    end
end
same_yaxes(AxesHandles_3)

%% Calculate single trace autocorrelation function

figure % plot 10 x 5 subplots. A_length = 10
num_examples = 5;
count3 = 0;
numLags = SAMPLES - 1;
lags = (1:numLags)';

% Force exponential fit to pass through autocorr = 1 at lag=0
x0 = 0;  % x coordinate of preestablished point
y0 = 1;  % y coordinate of preestablished point 
g = @(p,x)y0*exp(-p*(x-x0));
% f = fit(x,y,g)
t_decay = zeros(A_length, num_examples);
t_decay_ave = zeros(A_length, 1);

for i = 1:A_length
    this_X = X_csv{i}; % 1952 x ? matrix
    
    % Plot average trace in first column
    count3 = count3 + 1;
    subplot(A_length, num_examples + 2, count3)  
    y = mean(X_csv{i}, 2);
    plot(x, y)
    title(['Input A=', num2str(A_csv(i)), ' dB SPL'])
    xlabel('Time (ms)')
    ylabel('Voltage (nV)')
    
    % Plot ACR for averaged trace in second column
    count3 = count3 + 1;
    subplot(A_length, num_examples + 2, count3)  
    ta = acf(y, numLags, 0);
    bar(ta)
    yl=ylim;
    % Plot envelope
    [yupper,ylower] = envelope(ta);
%     f = fit(lags,yupper,'exp1');
    f = fit(lags, yupper, g, 'Lower', 0);
    t_decay_ave(i) = 1/f.p;
    hold on
%     plot(lags,yupper,lags,ylower,'linewidth',1.5)
    plot(lags,yupper,'linewidth',1.5)
    plot(f,lags,yupper)
    hLegend = findobj(gcf, 'Type', 'Legend');
    set(hLegend,'visible','off')
    hold off
    title('Average ACF')
    xlabel('Lag Length')
    ylabel('Autocorrelation')
    ylim(yl)
    
    % Loop over first <num_examples> single traces, and plot
    % autocorrelations in row
    for j = 1:num_examples
        y2 = this_X(:, j); 
        ta = acf(y2, numLags, 0);

        count3 = count3 + 1;
        subplot(A_length, num_examples + 2, count3)  

        % Plot ACF
        % Plot rejection region lines for test of individual autocorrelations
        % H_0: rho(tau) = 0 at alpha=.05
        bar(ta)
        yl=ylim;
        [yupper,ylower] = envelope(ta);
%         f = fit(lags,yupper,'exp1');
        f = fit(lags,yupper,g, 'Lower', 0);
        t_decay(i, j) = 1/f.p;
        hold on
%         plot(lags,yupper,lags,ylower,'linewidth',1.5)
        plot(lags,yupper,'linewidth',1.5)
        plot(f,lags,yupper)
        hLegend = findobj(gcf, 'Type', 'Legend');
        set(hLegend,'visible','off')
        hold off
        title(['T_{dec} = ', num2str(t_decay(i, j)), ' ms'])
        xlabel('Lag Length')
        ylabel('Autocorrelation')
        ylim(yl)
        
        p = numLags;
        N = max(size(y2));
        line([0 p+.5], (1.96)*(1/sqrt(N))*ones(1,2))
        line([0 p+.5], (-1.96)*(1/sqrt(N))*ones(1,2))

        % Some figure properties
        line_hi = (1.96)*(1/sqrt(N))+.05;
        line_lo = -(1.96)*(1/sqrt(N))-.05;
        bar_hi = max(ta)+.05 ;
        bar_lo = -max(ta)-.05 ;

        if (abs(line_hi) > abs(bar_hi)) % if rejection lines might not appear on graph
            axis([0 p+.60 line_lo line_hi])
        else
            axis([0 p+.60 bar_lo bar_hi])
        end
        set(gca,'YTick',[-1:.20:1])
        % set number of lag labels shown
        if (p<28 && p>4)
            set(gca,'XTick',floor(linspace(1,p,4)))
        elseif (p>=28)
            set(gca,'XTick',floor(linspace(1,p,8)))
        end
        set(gca,'TickLength',[0 0])
    end

end

%% Calculate area under curve of autocorrelation function
auc_td_ave = zeros(A_length, 1);
auc_td_single = cell(A_length,1); 
auc_td_single_ave = zeros(A_length, 1);
auc_td_single_std = zeros(A_length, 1);

t_decay_ave_rounded = round(t_decay_ave);

numLags = SAMPLES - 1;
for i = 1:A_length
    disp(['Working on ', num2str(A_csv(i)) ' of ', num2str(A_csv(A_length)), ' dB SPL...'])
    this_X = X_csv{i}; % 1952 x ? matrix
    y = mean(this_X, 2);
    ta = acf(y, numLags, 0);
    
    td = t_decay_ave_rounded(i);
    auc_td_ave(i) = trapz(abs(ta(1:td)));
    
    num_traces = size(this_X, 2);
    temp_auc_single = zeros(num_traces, 1);
    for j = 1:num_traces
        if mod(j, 50)==0
            disp(['  Single trace: ', num2str(j), ' of ', num2str(num_traces)])
        end
        y2 = this_X(:, j); 
        ta = acf(y2, numLags, 0);
        temp_auc_single(j) = trapz(abs(ta(1:td)));
    end
    
    % cache single trace AUC data
    auc_td_single{i} = temp_auc_single;
    auc_td_single_ave(i) = mean(temp_auc_single);
    auc_td_single_std(i) = std(temp_auc_single);
end

%% Calculate single trace power spectrum densities

figure % plot 10 x 5 subplots. A_length = 10
num_examples = 5;
count3 = 0;
numLags = SAMPLES - 1;
for i = 1:A_length
    this_X = X_csv{i}; % 1952 x ? matrix
    
    % Plot average trace in first column
    count3 = count3 + 1;
    subplot(A_length, num_examples + 2, count3)  
    y = mean(X_csv{i}, 2);
    plot(x, y)
    title(['Input A=', num2str(A_csv(i)), ' dB SPL'])
    xlabel('Time (ms)')
    ylabel('Voltage (nV)')
    
    % Plot ACR for averaged trace in second column
    count3 = count3 + 1;
    subplot(A_length, num_examples + 2, count3)  
    [pxx,f] = periodogram(y, [], [], SAMPLING_RATE);
    plot(f, pxx)
    title('Average PSD')
    xlim([0, 0.4*10^4])
    xlabel('Frequency')
    ylabel('PSD')
    
    % Loop over first <num_examples> single traces, and plot
    % autocorrelations in row
    for j = 1:num_examples
        y2 = this_X(:, j); 
        count3 = count3 + 1;
        subplot(A_length, num_examples + 2, count3)  

        % Plot PSD
        [pxx,f] = periodogram(y2, [], [], SAMPLING_RATE);
        plot(f, pxx)
        title(['Input A=', num2str(A_csv(i)), ' dB SPL'])
        xlabel('Frequency')
        ylabel('PSD')
        xlim([0, 0.4*10^4])
        
    end

end