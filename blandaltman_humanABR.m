% Bland altman plot of human ABR inter- and intra-subject variability of
% wave 1
%
% George Liu
% Last edit: 2/26/23
%
% Dependencies: none

PATH_ANALYSIS = 'D:\George-abr\human_abr\humanABR_analysis.xlsx';
FONT_SIZE = 14;
CIRCLE_SIZE = 36;
N_SUBJECTS = 5;

filename = PATH_ANALYSIS;
wave1amp_data = readmatrix(filename, 'Sheet', 'Wave 1 amplitude', 'Range', 'B2:C6');
wave1std_data = readmatrix(filename, 'Sheet', 'Wave 1 std', 'Range', 'B2:C6');
wave1lat_data = readmatrix(filename, 'Sheet', 'Wave 1 lat', 'Range', 'B2:C6');

%% Plot Bland-Altman of wave 1 amplitude
wave1amp_avg = nanmean(wave1amp_data, 2);
wave1amp_dif = abs(wave1amp_data(:,1) - wave1amp_data(:,2));

g1 = 1:N_SUBJECTS;
c = hsv(N_SUBJECTS);

% figure, scatter(wave1amp_avg, wave1amp_dif)
figure, gscatter(wave1amp_avg, wave1amp_dif, g1, c, '.', CIRCLE_SIZE)
axis equal
xlabel('Average wave 1 amplitude (nV)')
ylabel('Difference between two measurements (nV)')
set(gca,'FontSize', FONT_SIZE)

hold on
mean_dif_wave1 = nanmean(wave1amp_dif);
std_dif_wave1 = nanstd(wave1amp_dif);
yline(mean_dif_wave1, '-', sprintf('MEAN: %.0f', mean_dif_wave1), 'HandleVisibility','off')
yline(mean_dif_wave1 + std_dif_wave1, '--', sprintf('+SD: %.0f', mean_dif_wave1 + std_dif_wave1), 'HandleVisibility','off')
yline(mean_dif_wave1 - std_dif_wave1, '--', sprintf('-SD: %.0f', mean_dif_wave1 - std_dif_wave1), 'HandleVisibility','off')

%% Plot group scatter of wave 1 amplitude vs standard deviation
g2 = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5];
figure, gscatter(reshape(wave1amp_data.', 1, []), reshape(wave1std_data.', 1, []), g2, c, '.', CIRCLE_SIZE)
xlabel('Wave 1 amplitude (nV)')
ylabel('STD of single trace wave 1 amp (nV)')
set(gca,'FontSize', FONT_SIZE)

%% Plot Bland-Altman of wave 1 latency
wave1lat_avg = nanmean(wave1lat_data, 2);
wave1lat_dif = abs(wave1lat_data(:,1) - wave1lat_data(:,2));

figure, gscatter(wave1lat_avg, wave1lat_dif, g1, c, '.', CIRCLE_SIZE)
axis equal
xlabel('Average wave 1 latency (ms)')
ylabel('Difference between two measurements (ms)')
set(gca,'FontSize', FONT_SIZE)
