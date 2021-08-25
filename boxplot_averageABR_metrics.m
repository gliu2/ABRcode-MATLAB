% boxplot_averageABR_metrics
% script 
%
% 7/20/21 George Liu

%% Constants
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures

%% Select excel file without average ABR analysis results
[filename,path] = uigetfile('*.xlsx');
T = readtable(fullfile(path, filename));

% Add random jitter
rand_jitter = 300*rand(length(T.Frequency),1);
group_slide = [zeros(9,1); ones(9,1)]*1000;

% %% Scatter plots
% figure
% t = tiledlayout(2, 2);
% nexttile
% s = scatter(T.Frequency, T.Threshold, 'jitter','on', 'jitterAmount', 300);
% xlim([0, 35000])
% ylim([0, 90])
% xlabel('Frequency (Hz)')
% ylabel('Threshold (dB SPL)')
% title('ABR Threshold')
% s.LineWidth = 0.6;
% s.MarkerEdgeColor = 'b';
% s.MarkerFaceColor = [0 0.5 0.5];
% 
% gs = groupsummary(T,'Frequency','mean','Threshold');
% hold on
% plot(gs.Frequency, gs.mean_Threshold)
% hold off
% set(gca,'FontSize',10)
% 
% % Wave 1 amplitude
% nexttile
% s2 = scatter(T.Frequency, T.Amplitude, 'jitter','on', 'jitterAmount', 300);
% xlim([0, 35000])
% xlabel('Frequency (Hz)')
% ylabel('Amplitude (nV)')
% title('Wave 1 amplitude')
% s2.LineWidth = 0.6;
% s2.MarkerEdgeColor = 'b';
% s2.MarkerFaceColor = [0 0.5 0.5];
% 
% gs2 = groupsummary(T,'Frequency','mean','Amplitude');
% hold on
% plot(gs2.Frequency, gs2.mean_Amplitude)
% hold off
% set(gca,'FontSize',10)
% 
% % DPOAE threshold
% nexttile
% s3 = scatter(T.Frequency, T.DPOAEThreshold, 'jitter','on', 'jitterAmount', 300);
% xlim([0, 35000])
% ylim([0, 90])
% xlabel('Frequency (Hz)')
% ylabel('Threshold (dB)')
% title('DPOAE threshold')
% s3.LineWidth = 0.6;
% s3.MarkerEdgeColor = 'b';
% s3.MarkerFaceColor = [0 0.5 0.5];
% 
% gs3 = groupsummary(T,'Frequency','mean','DPOAEThreshold');
% hold on
% plot(gs3.Frequency, gs3.mean_DPOAEThreshold)
% hold off
% set(gca,'FontSize',10)
% 
% % Wave 1 latency
% nexttile
% s4 = scatter(T.Frequency, T.Latency, 'jitter','on', 'jitterAmount', 300);
% xlim([0, 35000])
% xlabel('Frequency (Hz)')
% ylabel('Latency (ms)')
% title('Wave 1 latency')
% s4.LineWidth = 0.6;
% s4.MarkerEdgeColor = 'b';
% s4.MarkerFaceColor = [0 0.5 0.5];
% 
% gs4 = groupsummary(T,'Frequency','mean','Latency');
% hold on
% plot(gs4.Frequency, gs4.mean_Latency)
% hold off
% set(gca,'FontSize',10)
% 
% 
% t.TileSpacing = 'tight';
% t.Padding = 'tight';
% 
% % Save figure
% [~, save_file, ~] = fileparts(filename);
% savefig_multiformat(gcf, SAVE_PATH, save_file)

%% Group scatter plots
figure
t = tiledlayout(2, 2);
nexttile
% gscatter(T.Frequency, T.Threshold, T.Calibrated);
gscatter(T.Frequency + rand_jitter + group_slide, T.Threshold, T.Calibrated);
% boxchart(T.Frequency, T.Threshold, 'GroupByColor', T.Calibrated)
xlim([0, 35000])
ylim([0, 90])
xlabel('Frequency (Hz)')
ylabel('Threshold (dB SPL)')
title('ABR Threshold')

gs = groupsummary(T,'Frequency','mean','Threshold');
hold on
plot(gs.Frequency, gs.mean_Threshold)
hold off
set(gca,'FontSize',10)
legend('Location','northeastoutside')

% Wave 1 amplitude
nexttile
% gscatter(T.Frequency, T.Amplitude, T.Calibrated);
gscatter(T.Frequency + rand_jitter + group_slide, T.Amplitude, T.Calibrated); 
% boxchart(T.Frequency, T.Amplitude, 'GroupByColor', T.Calibrated)
xlim([0, 35000])
xlabel('Frequency (Hz)')
ylabel('Amplitude (nV)')
title('Wave 1 amplitude')

gs2 = groupsummary(T,'Frequency','mean','Amplitude');
hold on
plot(gs2.Frequency, gs2.mean_Amplitude)
hold off
set(gca,'FontSize',10)
legend('Location','northeastoutside')

% DPOAE threshold
nexttile
% gscatter(T.Frequency, T.DPOAEThreshold, T.Calibrated);
gscatter(T.Frequency  + rand_jitter + group_slide, T.DPOAEThreshold, T.Calibrated);
% boxchart(T.Frequency, T.DPOAEThreshold, 'GroupByColor', T.Calibrated)
xlim([0, 35000])
ylim([0, 90])
xlabel('Frequency (Hz)')
ylabel('Threshold (dB)')
title('DPOAE threshold')

gs3 = groupsummary(T,'Frequency','mean','DPOAEThreshold');
hold on
plot(gs3.Frequency, gs3.mean_DPOAEThreshold)
hold off
set(gca,'FontSize',10)
legend('Location','northeastoutside')

% Wave 1 latency
nexttile
% gscatter(T.Frequency, T.Latency, T.Calibrated);
gscatter(T.Frequency + rand_jitter + group_slide, T.Latency, T.Calibrated);
% boxchart(T.Frequency, T.Latency, 'GroupByColor', T.Calibrated)
xlim([0, 35000])
xlabel('Frequency (Hz)')
ylabel('Latency (ms)')
title('Wave 1 latency')

gs4 = groupsummary(T,'Frequency','mean','Latency');
hold on
plot(gs4.Frequency, gs4.mean_Latency)
hold off
set(gca,'FontSize',10)
legend('Location','northeastoutside')

t.TileSpacing = 'tight';
t.Padding = 'tight';

% Save figure
[~, save_file, ~] = fileparts(filename);
savefig_multiformat(gcf, SAVE_PATH, save_file)

%%
figure, gscatter(T.Amplitude(1:9), T.Amplitude(10:end), T.Frequency(1:9))
hold on, plot(1:5000,1:5000), hold off
axis equal
xlim([0,6200])
ylim([0,6200])
xlabel('Wave 1 amplitude - 7/18/21 (nV)')
ylabel('Wave 1 amplitude - 7/14/21 (nV)')
title('Open ABR, uncalibrated')