%% jasaFigures_5_25_19.m
% 
% Create figures for JASA paper on inner product calculation of ABR
% single-trace ensemble threshold
%
% Run this script after "import_ABRcsv_folder.m", which loads Noor's
% single-trace ABR data.
%
% Dependencies: analyze_v2_ABR.m, same_yaxes.m, PTDetect.m 
% Last edit: 5/25/2019
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

%% 5-13-19: Plot averaged ABR waveforms on a relative scale -- "e.g. individual wave was normalized to the maximal wave within each ABR
% pattern", per Zhou et al. 2006 "Auditory brainstem responses in 10 inbred
% strains of mice" 

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
end
same_yaxes(AxesHandles_1)
    
%% 5-9-19: Plot averaged ABR RMS vs dB level (to see old threshold of 3*RMS no signal)
RMS_averaged_trace = sqrt(RMS2_averaged_trace);
figure('DefaultAxesFontSize', 20)
y_uV = RMS_averaged_trace*10^6; % in case want to change V -> uV y-axis units
NUM_STD_ABOVENOISE = 3;
thresh_1 = NUM_STD_ABOVENOISE*y_uV(1);
plot(A_csv, y_uV, '-o')
my_yline = yline(thresh_1, '--', ['THRESHOLD = ', num2str(NUM_STD_ABOVENOISE), '(RMS)_{0} = ', num2str(thresh_1), '\muV']);
title(['RMS of averaged ABR trace'], 'FontSize', 32)
xlabel('Stimulus level (dB SPL)', 'FontSize', 32)
ylabel('RMS (\muV)', 'FontSize', 32)
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


% Also plot 2*RMS_0 threshold line

%% 5-25-19: Threshold method 2: Find dB level when max peak-peak amplitude 
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
PEAK_THRESH = 0.5;
yave_cache = cell(A_length, 1);
yrel_cache = cell(A_length, 1);
peak_cache = cell(A_length, 1);
trough_cache = cell(A_length, 1);

figure
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
for i = 1:A_length
    y = yave_cache{i};
    P = peak_cache{i};
    T = trough_cache{i};
    
    % If peak first, then pair peak with trough and trough before
    if P(1) < T(1) % # troughs = # peaks OR # peaks-1
        peak2peak_pt = y(P(1:length(T))) - y(T); % peak then trough
        peak2peak_tp = y(P(2:end)) - y(T(1:length(P)-1)); % trough then peak
    else % T(1) < P(1) % # troughs = # peaks OR # peaks+1
        peak2peak_pt = y(P(1:length(T)-1)) - y(T(2:end));
        peak2peak_tp = y(P) - y(T(1:length(P)));
    end
    
    % Find maximum peak-to-peak amplitude
    max_p2p(i) = max([peak2peak_pt; peak2peak_tp]);
end

%% 5-25-19: Plot maximum peak-to-peak amplitude vs dB level
figure('DefaultAxesFontSize', 20)
RMS_averaged_trace_uV = RMS_averaged_trace*10^6; % change standard deviation (RMS) units V -> uV y-axis units
y_uV = max_p2p*10^6; % in case want to change V -> uV y-axis units
NUM_STD_ABOVENOISE = 4; % 4 or 5
thresh_up = NUM_STD_ABOVENOISE*RMS_averaged_trace_uV(1);
thresh_1 = thresh_up + y_uV(1);
plot(A_csv, y_uV, '-o')
my_yline = yline(thresh_1, '--', ['THRESHOLD = ', num2str(NUM_STD_ABOVENOISE), '(RMS)_{0} + V_0^{p-p} = ', num2str(thresh_up), ' \muV + ', num2str(y_uV(1)), ' \muV = ', num2str(thresh_1), ' \muV']);
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