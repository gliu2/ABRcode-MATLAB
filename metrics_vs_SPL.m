%% metrics_vs_SPL.m
% 
% Plot various metrics of "signal" in ABR single traces to see if any
% correlate with input tone pip amplitude.
%
% Run this script after "import_ABRcsv_folder.m", which loads Noor's
% single-trace ABR data.
%
% Assesses: <V(tn)> - averaged ABR voltage across all single-traces at the
% same time-point in each single trace.
%
% Dependencies: none
% Last edit: 3/18/2019
%
% Author: George Liu

SAMPLES = size(X_csv{i}, 1);
SAMPLING_RATE = 200000; % Hz
dt = 1/SAMPLING_RATE * 1000; % ms
T_total = (SAMPLES-1)*dt; % ms
x = 0:dt:T_total;
% y1 = sin(x);
% y2 = cos(x);
% y3= tan(x);
% y4=1./cos(x);

%% Plot averaged ABR waveforms -- the current standard for visually determining ABR threshold

RMS2_averaged_trace = zeros(A_length, 1);
X_csv_avg_all = zeros(SAMPLES, A_length);
sigma_v_tn = zeros(SAMPLES, A_length); % SAMPLES x # SPL's matrix of amplitude standard deviation across single traces at fixed time point
RMS_I = zeros(A_length, 1);
sigma_RMSI = zeros(A_length, 1);

figure
for i = 1:A_length
    y = mean(X_csv{i}, 2);
    
    subplot(ceil(A_length/2),2,i)  
    plot(x, y)
    title(['Input A=', num2str(A_csv(i)), ' dB SPL'])
    xlabel('Time (ms)')
    ylabel('Voltage (nV)')
    
    % Append RMS squared
    RMS2_averaged_trace(i, 1) = analyze_v2_ABR(y);
    
    % Calculate std at fixed time point across single traces
    std_this_v = std(X_csv{i}, 0, 2); % usage: S = std(A,w,dim) 
    sigma_v_tn(:, i) = std_this_v;
    
    % Calculate mean RMS of single traces for each intensity stimulus
    RMS_this_v = sqrt(analyze_v2_ABR(X_csv{i})); % m x 1 vector
    RMS_I(i) = mean(RMS_this_v);
    sigma_RMSI(i) = std(RMS_this_v);
    
    % Cache average ABR plot
    X_csv_avg_all(:, i) = y;
end

%% Plot <V(tn)> - averaged ABR voltage across all single-traces at the
% same time-point in each single trace
figure
xx = repmat(x',1,A_length);
yy = repmat(A_csv', SAMPLES, 1);

h1 = surf(xx, yy, X_csv_avg_all*10^6);
set(h1,'LineStyle','none')
xlabel('Time (ms)')
ylabel('dB SPL')
zlabel('\muV')
title('<V(t_n)>')
colorbar

%% Sigma_v(tn)
figure
h2 = surf(xx, yy, sigma_v_tn*10^6);
set(h2,'LineStyle','none')
xlabel('Time (ms)')
ylabel('dB SPL')
zlabel('\muV')
title('\sigma_v(tn)')
colorbar

%% Coefficient of variation <V(tn)>/sigma(tn)
figure
h3 = surf(xx, yy, X_csv_avg_all./sigma_v_tn);
set(h3,'LineStyle','none')
xlabel('Time (ms)')
ylabel('dB SPL')
zlabel('a.u.')
title('<V(tn)>/\sigma_v(tn)')
colorbar

%% <RMS_I>
figure
plot(A_csv, RMS_I*10^6)
xlabel('dB SPL')
ylabel('<RMS_I> (\muV)')
title('Single trace <RMS>')

%% sigma_RMS_I
figure
plot(A_csv, sigma_RMSI*10^6)
xlabel('dB SPL')
ylabel('\sigma_{RMS_I} (\muV)')
title('Single trace sigma_RMS')

%% <RMS_I>/sigma_RMS_I
figure
plot(A_csv, RMS_I./sigma_RMSI)
xlabel('dB SPL')
ylabel('<RMS_I>/\sigma_{RMS_I} (a.u.)')
title('Single trace CV^{-1}')
