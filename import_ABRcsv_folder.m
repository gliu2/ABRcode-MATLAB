% Read ABR responses for all CSV spreadsheets in folder
% Dependencies: natsortorder.m

%% Get folder path
selpath = uigetdir();
dir_csv = dir(fullfile(selpath, '*.csv'));

%% Import data
A_length = length(dir_csv);
X_csv = cell(A_length, 1);
A_csv = zeros(A_length, 1);
freq_csv = zeros(A_length, 1);

%Obtain natural order of file names
filenames = cell(A_length, 1);
for i = 1:A_length
    filenames{i} = dir_csv(i).name;
end
filename_natorder = natsortfiles(filenames);

for i = 1:A_length
%     filename = dir_csv(i).name;
%     [X_csv{i}, A_csv(i), freq_csv(i)] = import_ABRcsv(filename, selpath);
    [X_csv{i}, A_csv(i), freq_csv(i)] = import_ABRcsv(filename_natorder{i}, selpath);
end

%% Get new and old psychometric curves
X_noise_csv = X_csv{1};
X_noise_avg_csv = mean(X_noise_csv, 2); % SAMPLES x 1 vector
std_old = std(X_noise_avg_csv);

MULT = 2; % # std of noise
RMS_old_cutoff = MULT*std_old;
Vd_old = RMS_old_cutoff.^2; % old detection threshold for RMS squared

m = zeros(size(A_csv));
Phit_old = zeros(size(A_csv));
Phit_new = zeros(size(A_csv));
for i = 1:A_length
    % Old method
    X_csv_avg = mean(X_csv{i}, 2); % SAMPLES x 1 vector
    RMS2_old = analyze_v2_ABR(X_csv_avg);
    Phit_old(i) = sum(RMS2_old > Vd_old);
    
    % New method
    m(i) = size(X_csv{i}, 2);
    RMS_new_cutoff = RMS_old_cutoff * sqrt(m(i));
    Vd_new = RMS_new_cutoff.^2; % new detection threshold for RMS squared
    RMS2_new = analyze_v2_ABR(X_csv{i});
    Phit_new(i) = sum(RMS2_new > Vd_new)/m(i);
end

% %% Get experimental ABR thresholds
% A_exp_old = ZeroGL(A_csv, Phit_old-0.5);  % not sure why also returns  9.7000
% if isempty(A_exp_old)
%     A_exp_old = A(end); % if no hits, set threshold to maximum possible value
% else 
%     A_exp_old = A_exp_old(1);
% end
% A_exp_new = ZeroGL(A_csv, Phit_new-0.5);  % not sure why also returns  9.7000
% if isempty(A_exp_new)
%     A_exp_new = A(end); % if no hits, set threshold to maximum possible value
% else 
%     A_exp_new = A_exp_new(1);
% end
% all_thresholds = [A_exp_old, A_exp_new]; % for plotting

%% Plot psychometric-like curve
h=figure; plot(A_csv, Phit_new, 'LineWidth', 1, 'Color', 'red')
hold on
plot(A_csv, Phit_old, 'LineWidth', 1, 'Color', 'blue')
% plot(A_theoretical, hitrate_theoretical, 'LineWidth', 1, 'Color', 'green')
% plot(all_thresholds, ones(size(all_thresholds))/2, 'pg') % thresholds
hold off
xlabel('SPL_{TH} (dB)')
% xlabel('Signal ABR amplitude (\muV)', 'FontSize', 32)
ylabel('Probability of hit', 'FontSize', 32)
% title(['Psychometric curve for signal ABR detection (Vd=', num2str(Vd), ' {\muV}^2, FP=', num2str(FPR), ')'])
title(['Psychometric curve (Vd=', num2str(Vd_old), ' {nV}^2)'], 'FontSize', 32)
box off % remove ticks on top and right borders of plot
% axis tight % makes edges of data flush with left and right borders of plot
legend(['New method (m=', num2str(mean(m)), ' single traces)'],'Old method (1 averaged trace)', ...
    'location','northeast', 'FontSize', 16)
legend boxoff % remove box around legend
set(gca,'FontSize',20)