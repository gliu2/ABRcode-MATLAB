% plot_synapse_counts.m
% Quick script to plot synapse count data (manually stored to xy variable)
% Run plot_ABRthreshold_all_scatter.m to initialize constants first
%
% George Liu
% Last edit 12/20/21
% Created 11/21/21 

% Initialize constants
P_CRIT = 0.05;
LINEWIDTH = 2;
FONTSIZE = 14;
SYNAPSE_FILEPATH = 'd:\users\admin\Documents\George\synapse_counts_11-21-21.xlsx';
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures
USE_HIGHQUALITY = true; % hard code change

%% Import synapse count data from Excel
opts_isnoise = spreadsheetImportOptions('Sheet', 2, 'DataRange', 'A2:A167', 'VariableTypes', 'double');
opts_noise = spreadsheetImportOptions('Sheet', 2, 'DataRange', 'B2:B167', 'VariableTypes', 'double');
opts_frequency = spreadsheetImportOptions('Sheet', 2, 'DataRange', 'C2:C167', 'VariableTypes', 'double');
% opts_synapse = spreadsheetImportOptions('Sheet', 2, 'DataRange', 'Q2:Q167', 'VariableTypes', 'double');
opts_synapse = spreadsheetImportOptions('Sheet', 2, 'DataRange', 'R2:R167', 'VariableTypes', 'double'); % ignore synapse counts with less than 120 colocalized puncta, as these values are too low
opts_ohc = spreadsheetImportOptions('Sheet', 2, 'DataRange', 'S2:S167', 'VariableTypes', 'double'); 
opts_ohc_good = spreadsheetImportOptions('Sheet', 2, 'DataRange', 'T2:T167', 'VariableTypes', 'double'); % don't use low quality OHC count data, as determined from notes in excel sheet
opts_quality_ctbp2 = spreadsheetImportOptions('Sheet', 2, 'DataRange', 'W2:W167', 'VariableTypes', 'double'); % quality of puncta labeling in image, ramge 0-3, more stars is better
opts_quality_homer = spreadsheetImportOptions('Sheet', 2, 'DataRange', 'X2:X167', 'VariableTypes', 'double'); % quality of puncta labeling in image, ramge 0-3, more stars is better

is_noise = readmatrix(SYNAPSE_FILEPATH, opts_isnoise);
is_noise = logical(is_noise);
noise_level = readmatrix(SYNAPSE_FILEPATH, opts_noise);
freq = readmatrix(SYNAPSE_FILEPATH, opts_frequency);
coloc_ihc_ratio = readmatrix(SYNAPSE_FILEPATH, opts_synapse);
ohc_ratio = readmatrix(SYNAPSE_FILEPATH, opts_ohc);
ohc_ratio_good = readmatrix(SYNAPSE_FILEPATH, opts_ohc_good);
quality_ctbp2 = readmatrix(SYNAPSE_FILEPATH, opts_quality_ctbp2);
quality_homer = readmatrix(SYNAPSE_FILEPATH, opts_quality_homer);

%% 4-18-22 Import synapse count file name metadata
IND_COL_90DB = 24; % column number of 90 dB amplitude values

opts_filename = spreadsheetImportOptions('Sheet', 2, 'DataRange', 'D2:D167', 'VariableTypes', 'string');
filenames = readmatrix(SYNAPSE_FILEPATH, opts_filename);

% For each synapse count, import abr wave 1 amplitude
% all_date = finalTable_metadata.('Date');
% all_date1 = finalTable_metadata.('Date_1');
% all_date2 = finalTable_metadata.('Date_2');
% all_date3 = finalTable_metadata.('Date_3');
metric_selector = strcmpi(finalTable_metadata.('Metric'), 'Wave1amp');

% Obtain the wave 1 amplitude pre, post24h and post1w exposure for each
% synapse count
n_files = numel(filenames);

wave1amp_pre = NaN(n_files, 1);
wave1amp_post24h = NaN(n_files, 1);
wave1amp_post1w = NaN(n_files, 1);
for i=1:n_files
    filename_i = filenames(i);
    [date, name, studytype, side, metadata] = get_mousefile_metadata(filename_i);
    
    % "studytype" 3rd input is actually side - rename variables
    section = side;
    side = studytype;
    
    is_side = strcmpi(finalTable_metadata.('Side'), side);
    is_name = strcmpi(finalTable_metadata.('Mouse_name'), name);
    is_freq = finalTable_metadata.('Frequency')== freq(i)*1000;
%     is_date = strcmpi(finalTable_metadata.('Date'), date);

    all_filtered = is_side & is_name & is_freq & metric_selector;
    
    this_data = finalTable_metadata(all_filtered, :);
    
    % Get wave 1 amplitude values before, 24h, and 1w after noise exposure
    date1 = NaN;
    date2 = NaN;
    date3 = NaN;
    if numel(this_data.Date_1) > 0
        date1 = this_data.Date_1(1);
    end
    if numel(this_data.Date_2) > 0
        date2 = this_data.Date_2(1);
    end
    if numel(this_data.Date_3) > 0
        date3 = this_data.Date_3(1);
    end
   
        
    % Ensure dates are not cell arrays
    if iscell(date1)
        date1 = date1{1};
    end
    if iscell(date2)
        date2 = date2{1};
    end
    if iscell(date3)
        date3 = date3{1};
    end
    
    % Ensure dates are numbers
    if isstring(date1) || ischar(date1)
        date1 = str2double(date1);
    end
    if isstring(date2) || ischar(date2)
        date2 = str2double(date2);
    end
    if isstring(date3) || ischar(date3)
        date3 = str2double(date3);
    end


    % Iterate over every ABR date for that mouse ear
    for j = 1:size(this_data, 1)
        % Date of ABR
        this_date = this_data.Date(j);
        
        if isstring(this_date)
            this_date = str2double(this_date);
        end
        
        % Get wave 1 amplitude for this date
        this_wave1amp = this_data{j, IND_COL_90DB};
        
        % Determine date of wave 1 amplitude relative to noise exposure
        if this_date == date1
            wave1amp_pre(i) = this_wave1amp;
        elseif this_date == date2
            wave1amp_post24h(i) = this_wave1amp;
        elseif this_date == date3 
            wave1amp_post1w(i) = this_wave1amp;
        end
    end
end

%% 4-21-22 Plot wave 1 amplitude vs synapse count

all_x = [coloc_ihc_ratio(~is_noise); coloc_ihc_ratio(is_noise)];
all_y1 = [wave1amp_post1w(~is_noise); wave1amp_post1w(is_noise)];
all_y2 = [wave1amp_post1w(~is_noise) - wave1amp_pre(~is_noise); wave1amp_post1w(is_noise) - wave1amp_pre(is_noise)];

xlim_this = [min(all_x), max(all_x)];
ylim_this1 = [min(all_y1), max(all_y1)];
ylim_this2 = [min(all_y2), max(all_y2)];

figure
t=tiledlayout(2,2);

% Control mice
nexttile
gscatter(coloc_ihc_ratio(~is_noise), wave1amp_post1w(~is_noise), noise_level(~is_noise))
title('Control')
xlabel('# Colocalized synapses per IHC')
ylabel('Wave 1 amplitude (nV) at 1 week')
xlim(xlim_this)
ylim(ylim_this1)

% Noise exposed mice
nexttile
gscatter(coloc_ihc_ratio(is_noise), wave1amp_post1w(is_noise), noise_level(is_noise))
title('Noise exposed')
xlabel('# Colocalized synapses per IHC')
ylabel('Wave 1 amplitude (nV) at 1 week')
xlim(xlim_this)
ylim(ylim_this1)

% Control mice
nexttile
gscatter(coloc_ihc_ratio(~is_noise), wave1amp_post1w(~is_noise) - wave1amp_pre(~is_noise), noise_level(~is_noise))
title('Control')
xlabel('# Colocalized synapses per IHC')
ylabel('\Delta Wave 1 amplitude, 1 week - baseline (nV)')
xlim(xlim_this)
ylim(ylim_this2)

% Noise exposed mice
nexttile
gscatter(coloc_ihc_ratio(is_noise), wave1amp_post1w(is_noise) - wave1amp_pre(is_noise), noise_level(is_noise))
title('Noise exposed')
xlabel('# Colocalized synapses per IHC')
ylabel('\Delta Wave 1 amplitude, 1 week - baseline (nV)')
xlim(xlim_this)
ylim(ylim_this2)

t.Padding = 'compact';
t.TileSpacing = 'compact';

%% Plot synapse counts (ratio of synapses to IHCs) per frequency and noise exposure level
n_noise = sum(is_noise);
n_control = length(is_noise) - n_noise;
freqs_unique = unique(freq);
noise_levels_unique = unique(noise_level);
n_freq = numel(freqs_unique);
n_noise_levels = numel(noise_levels_unique);

ishighquality = quality_ctbp2>1 & quality_homer>1; % filter by quality of image data if desired

% group_labels = {['Control, n=', num2str(n_control)], ['96 dB, n=', num2str(n_noise)]};
% figure, boxplot(coloc_ihc_ratio, [is_noise, freq, noise_level], 'Labels', group_labels)

% figure, boxplot(coloc_ihc_ratio, [noise_level, freq, is_noise], 'Labels', {'Noise level (dB)', 'Freq (kHz)', 'Is noise'})
groups = [noise_level, freq, is_noise];
fig = figure; 
if USE_HIGHQUALITY
    boxplot(coloc_ihc_ratio(ishighquality), groups(ishighquality, :), 'Notch','off')
else
    boxplot(coloc_ihc_ratio, groups, 'Notch','off')
end
title('Synapse counts')
ylabel('# colocalized synapses per IHC')
set(gca,'FontSize', FONTSIZE) % Set text size of current axis

% hold on
% s=scatter(is_noise + 1.2, coloc_ihc_ratio, 'filled');
% s.AlphaData = 0.1*ones(size(coloc_ihc_ratio)); % set transparency of filled circle marks
% s.MarkerFaceAlpha = 'flat';
% 
% [~, p_synapse] = ttest2(coloc_ihc_ratio(~is_noise), coloc_ihc_ratio(logical(is_noise)), 'Tail', 'right');
% if p_synapse < 0.05
%     text(2.3, mean(coloc_ihc_ratio(logical(is_noise))), '*', 'FontSize', FONTSIZE, 'FontWeight', 'bold')
% end
% hold off
% 

% % Maximize figure window size
% fig.WindowState = 'maximized';

% Save figure to results folder
save_filename = ['synapse_boxplot_usehighqualdata_', num2str(USE_HIGHQUALITY)];
savefig_multiformat(gcf, SAVE_PATH, save_filename)

%% Plot synapse count vs frequency, pooling all controls together
n_groups = n_noise_levels + 1;
all_noise_levels = vertcat(0, noise_levels_unique);
synapse_mean = zeros(n_groups, n_freq);
synapse_ste = zeros(n_groups, n_freq);
n_counts = zeros(n_groups, n_freq);
synapse_counts = cell(n_groups, n_freq);

% For control and each noise exposure level, plot separate line
colors = ['k', 'b', 'g', 'r'];
fig = figure;
for i=1:n_groups
    this_noise_level = all_noise_levels(i);
    
    % Iterate over frequencies
    for j=1:n_freq
        this_freq = freqs_unique(j);
        
        if USE_HIGHQUALITY % use only synapse data points from images with >1 stars quality for ctbp2 and homer staining
            if this_noise_level==0 % control
                this_synapse_counts = coloc_ihc_ratio(~is_noise & freq==this_freq & ishighquality);
            else % noise
                this_synapse_counts = coloc_ihc_ratio(is_noise & noise_level==this_noise_level & freq==this_freq & ishighquality);
            end
        else % use all data points regardliness of image quality
            if this_noise_level==0 % control
                this_synapse_counts = coloc_ihc_ratio(~is_noise & freq==this_freq);
            else % noise
                this_synapse_counts = coloc_ihc_ratio(is_noise & noise_level==this_noise_level & freq==this_freq);
            end
        end

        n_counts(i, j) = sum(~isnan(this_synapse_counts));
        synapse_mean(i, j) = nanmean(this_synapse_counts);
        synapse_ste(i,j) = nanstd(this_synapse_counts)/sqrt(n_counts(i, j));
        synapse_counts{i,j} = this_synapse_counts;
    end
    
    % Plot
    errorbar(freqs_unique, synapse_mean(i, :), synapse_ste(i, :), colors(i), 'LineWidth', LINEWIDTH)
    hold on
end

% Plot asterisk above plots if noise synapse count is significantly less
% than control
for i=2:n_groups
    for j=1:n_freq
        [~, p] = ttest2(synapse_counts{1,j}, synapse_counts{i,j}, 'Tail', 'right'); % compare control and noise synapse means
        if p < P_CRIT
            text(freqs_unique(j), synapse_mean(i, j), '*', 'Color', colors(i), 'FontSize', FONTSIZE*2, 'FontWeight', 'bold')
        end
    end
end
hold off
legend_labels = {'Control', '96 dB', '97 dB', '99 dB'};
for k=1:n_groups
        legend_labels{k} = [legend_labels{k}, ', n=', num2str(n_counts(k, 1)), '/', num2str(n_counts(k, 2)), '/', num2str(n_counts(k, 3))];
end
legend(legend_labels, 'Location', 'southwest')
xlabel('Frequency (kHz)')
ylabel('# colocalized synapses per IHC')
ylim([0, 23])
title('Synapse counts');
set(gca,'FontSize', FONTSIZE) % Set text size of current axis
% % Maximize figure window size
% fig.WindowState = 'maximized';
            
% Save figure to results folder
save_filename = ['synapse_vs_freq_usehighqualdata_', num2str(USE_HIGHQUALITY)];
savefig_multiformat(gcf, SAVE_PATH, save_filename)


%% 12-28-21: Plot synapse count vs frequency, for 96 and 97 dB groups only. 2 plots. Do not pool controls.
TWO_NOISE_LEVELS = [96, 97]; % dB
n_noise_levels = numel(TWO_NOISE_LEVELS);

n_groups = n_noise_levels*2;
synapse_mean = zeros(n_groups, n_freq);
synapse_ste = zeros(n_groups, n_freq);
n_counts = zeros(n_groups, n_freq);
synapse_counts = cell(n_groups, n_freq);

% For control and each noise exposure level, plot separate line
colors = ['k', 'r'];
linestyles = {'-', '--'};
fig = figure;
count = 0;
for i=1:n_noise_levels
    this_noise_level = TWO_NOISE_LEVELS(i);
    
    for iscontrol=[true, false]
        count = count + 1;
        
        % Iterate over frequencies
        for j=1:n_freq
            this_freq = freqs_unique(j);

            if USE_HIGHQUALITY % use only synapse data points from images with >1 stars quality for ctbp2 and homer staining
                if iscontrol % control
                    this_synapse_counts = coloc_ihc_ratio(~is_noise & noise_level==this_noise_level & freq==this_freq & ishighquality);
                else % noise
                    this_synapse_counts = coloc_ihc_ratio(is_noise & noise_level==this_noise_level & freq==this_freq & ishighquality);
                end
            else % use all data points regardliness of image quality
                if iscontrol % control
                    this_synapse_counts = coloc_ihc_ratio(~is_noise & noise_level==this_noise_level & freq==this_freq);
                else % noise
                    this_synapse_counts = coloc_ihc_ratio(is_noise & noise_level==this_noise_level & freq==this_freq);
                end
            end

            n_counts(count, j) = sum(~isnan(this_synapse_counts));
            synapse_mean(count, j) = nanmean(this_synapse_counts);
            synapse_ste(count,j) = nanstd(this_synapse_counts)/sqrt(n_counts(count, j));
            synapse_counts{count,j} = this_synapse_counts;
        end
        
        % Plot
        if iscontrol
            this_linestyle = linestyles{1};
        else
            this_linestyle = linestyles{2};
        end
        errorbar(freqs_unique, synapse_mean(count, :), synapse_ste(count, :), colors(i), 'LineStyle', this_linestyle, 'LineWidth', LINEWIDTH)
        hold on
        
    end
    
end

% Plot asterisk above plots if noise synapse count is significantly less
% than control
for i=1:n_noise_levels
    for j=1:n_freq
        [~, p] = ttest2(synapse_counts{2*i-1,j}, synapse_counts{2*i,j}, 'Tail', 'both'); % compare control and noise synapse means
        if p < P_CRIT
            text(freqs_unique(j), synapse_mean(2*i, j), '*', 'Color', colors(i), 'FontSize', FONTSIZE*2, 'FontWeight', 'bold')
        end
    end
end
hold off
legend_labels = {'96 dB, control', '96 dB, noise', '97 dB, control', '97 dB, noise'};
for k=1:n_noise_levels*2
        legend_labels{k} = [legend_labels{k}, ', n=', num2str(n_counts(k, 1)), '/', num2str(n_counts(k, 2)), '/', num2str(n_counts(k, 3))];
end
legend(legend_labels, 'Location', 'southwest')
xlabel('Frequency (kHz)')
ylabel('# colocalized synapses per IHC')
ylim([0, 23])
title('Synapse counts');
set(gca,'FontSize', FONTSIZE) % Set text size of current axis
% % Maximize figure window size
% fig.WindowState = 'maximized';
            
% Save figure to results folder
save_filename = ['synapse_vs_freq_unpooled_usehighqualdata_', num2str(USE_HIGHQUALITY)];
savefig_multiformat(gcf, SAVE_PATH, save_filename)

%% Plot OHC : IHC ratio vs frequency
n_groups = n_noise_levels + 1;
all_noise_levels = vertcat(0, noise_levels_unique);
ohc_mean = zeros(n_groups, n_freq);
ohc_ste = zeros(n_groups, n_freq);
n_counts = zeros(n_groups, n_freq);
ohc_counts = cell(n_groups, n_freq);

% For control and each noise exposure level, plot separate line
colors = ['k', 'b', 'g', 'r'];
fig = figure;
for i=1:n_groups
    this_noise_level = all_noise_levels(i);
    
    % Iterate over frequencies
    for j=1:n_freq
        this_freq = freqs_unique(j);
        
        if USE_HIGHQUALITY % use only synapse data points from images with >1 stars quality for ctbp2 and homer staining
            if this_noise_level==0 % control
                this_ohc_counts = ohc_ratio_good(~is_noise & freq==this_freq);
            else % noise
                this_ohc_counts = ohc_ratio_good(is_noise & noise_level==this_noise_level & freq==this_freq);
            end
        else % use all data points regardliness of image quality
            if this_noise_level==0 % control
                this_ohc_counts = ohc_ratio(~is_noise & freq==this_freq);
            else % noise
                this_ohc_counts = ohc_ratio(is_noise & noise_level==this_noise_level & freq==this_freq);
            end
        end

        n_counts(i, j) = sum(~isnan(this_ohc_counts));
        ohc_mean(i, j) = nanmean(this_ohc_counts);
        ohc_ste(i,j) = nanstd(this_ohc_counts)/sqrt(n_counts(i, j));
        ohc_counts{i,j} = this_ohc_counts;
    end
    
    % Plot
    errorbar(freqs_unique, ohc_mean(i, :), ohc_ste(i, :), colors(i), 'LineWidth', LINEWIDTH)
    hold on
end

% Plot asterisk above plots if noise synapse count is significantly less
% than control
for i=2:n_groups
    for j=1:n_freq
        [~, p] = ttest2(ohc_counts{1,j}, ohc_counts{i,j}, 'Tail', 'both'); % compare control and noise synapse means
        if p < P_CRIT
            text(freqs_unique(j), ohc_mean(i, j), '*', 'Color', colors(i), 'FontSize', FONTSIZE*2, 'FontWeight', 'bold')
        end
    end
end
hold off
legend_labels = {'Control', '96 dB', '97 dB', '99 dB'};
for k=1:n_groups
        legend_labels{k} = [legend_labels{k}, ', n=', num2str(n_counts(k, 1)), '/', num2str(n_counts(k, 2)), '/', num2str(n_counts(k, 3))];
end
legend(legend_labels, 'Location', 'southwest')
xlabel('Frequency (kHz)')
ylabel('# OHCs per IHC')
ylim([1, 4])
title('OHC counts');
set(gca,'FontSize', FONTSIZE) % Set text size of current axis
% % Maximize figure window size
% fig.WindowState = 'maximized';
            
% Save figure to results folder
save_filename = ['ohc_vs_freq_usehighqualdata_', num2str(USE_HIGHQUALITY)];
savefig_multiformat(gcf, SAVE_PATH, save_filename)
