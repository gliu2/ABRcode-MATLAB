function [thresh, metric] = get_thresh_CAPABR_D_allfreq(X_csv_cache, A_descending, all_freqs, varargin)
% Calculate threshold for all frequencies using a metric of D', i.e. Cohen's D. 
% Uses single trial inner product. 
%
% Input: X_csv_cache - sorted cell array of single trial ABR traces. Cell of size (num levels, n_freq), 
%                   Each cell entry is a matrix of size (SAMPLES, m_traces), for set of single ABR
%                   traces at same dB level. Columns are individual trace data.
%        A_descending - vector of stimulus levels (dB) in descending order
%        all_freqs - column vector of frequencies corresponding to columns
%                   of X_csv_cache
%        PLOT_FIGURES (optional) - boolean, specify whether to plot figures
%
% output: 
%       thresh - threshold level (dB) calculated using metric. Matrix of
%                 size (n_freq, 1).
%       metric - D' compared with 0 dB response. Matrix of size
%               (num_levels, n_freq).
%
% Last edit: 8/25/23 George Liu
%
% Dependencies: get_thresh_CAPABR_D.m

%% Constant
cutoff = 1;
PLOT_FIGURES = true;

if nargin==4
    PLOT_FIGURES = varargin{1};
end

num_levels = length(A_descending);
n_freq = length(all_freqs);
thresh = zeros(n_freq, 1);
metric = zeros(num_levels, n_freq);
for i = 1:n_freq 
    [this_thresh, this_metric] = get_thresh_CAPABR_D(X_csv_cache(:, i), A_descending, cutoff, false);
    
    thresh(i) = this_thresh;
    metric(:, i) = this_metric;
end

%% Plot threshold
if ~PLOT_FIGURES
    return
end

%Plot average versus SPL
figure('DefaultAxesFontSize', 20)
ylim_bounds = [min([min(metric, [], "all"), 0]), 1.05*max(metric, [], "all")];
hold on
for i = 1:n_freq
    plot(A_descending, metric(:, i), 'LineWidth', 3);
    
%     % Draw vertical line at threshold amplitude 
%     this_thresh = thresh(i);
%     if ~isnan(this_thresh)
%         line([this_thresh this_thresh], [ylim_bounds(1), interp1(A_descending, metric(:, i), this_thresh, 'linear')], 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--'); % vertical line at exact threshold
%     else
%         disp('Warning: No ABR threshold detected!')
%     end
end

yline_label = ['Cutoff: D = ', num2str(cutoff)];
% my_yline = yline(cutoff, '--', ['Cutoff @ D = ', num2str(cutoff)]);
my_yline = yline(cutoff, '--');

all_freqs_khz = all_freqs/1000;
freq_units = cell(n_freq, 1);
freq_units(:) = {'kHz:'};
dB = cell(n_freq, 1);
dB(:) = {'dB'};

freq_cell = cellstr(num2str(all_freqs_khz));
thresh_cell = cellstr(num2str(round(thresh, 1)));
legend_labels = cellfun(@(w, x, y, z) [w ' ' x ' ' y ' ' z], freq_cell, freq_units, thresh_cell, dB, 'un', 0);
legend_labels = [legend_labels; {yline_label}];
legend(legend_labels, 'Location', 'northwest')

ylabel('Cohen D (a.u.)', 'FontSize', 24)
title('Threshold all frequencies', 'FontSize', 24)
xlabel('Stimulus level (dB SPL)', 'FontSize', 24)
set(my_yline, 'FontSize', 18)
set(gca,'TickDir','out');
set(gca,'box','off')
xlim([A_descending(end)-5, A_descending(1)+5]);
ylim(ylim_bounds);

end
