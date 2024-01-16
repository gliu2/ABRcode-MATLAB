function [thresh, metric, dist_innerprod] = get_thresh_CAPABR_D(M_sorted, A_descending, varargin)
% Calculate threshold using a metric of D', i.e. Cohen's D. 
% Uses single trial inner product. 
%
% Input: M_sorted - sorted cell array of single trial ABR traces. Cell of size (num levels, 1), 
%                   Each cell entry is a matrix of size (SAMPLES, m_traces), for set of single ABR
%            traces at same dB level. Columns are individual trace data.
%        A_descending - vector of stimulus levels (dB) in descending order
%        cutoff (optional) - cutoff value of D' for determining threshold.
%        PLOT_FIGURES (optional) - boolean, specify whether to plot figures
%
% output: 
%       thresh - threshold level (dB) calculated using metric
%       metric - D' compared with 0 dB response
%
% Last edit: 8/30/23 George Liu
%
% Dependencies: cohen_d_unpooled.m

%% Constant
cutoff = 1;
MONOTONIC_INCREASE_BELOW_THIS_D = 5;
PLOT_FIGURES = true;

if nargin==3
    cutoff = varargin{1};
elseif nargin == 4
    cutoff = varargin{1};
    PLOT_FIGURES = varargin{2};
end

%% Calculate ABR threshold and level function using max peak-peak amplitude
%Find maximum peak-to-peak amplitude in average ABR at each dB level
A_length = length(A_descending);
metric = zeros(A_length, 1);
dist_innerprod = zeros(A_length, size(M_sorted{1}, 2));

% Calculate distribution of inner product at 0 dB
dist_innerprod(A_length, :) = get_innerprod_matrix_ABR(M_sorted{A_length});

for i=1:A_length 
    dist_innerprod(i, :) = get_innerprod_matrix_ABR(M_sorted{i});
    
    % Calculate Cohen's D
    mean1 = mean(dist_innerprod(i, :));
    var1 = var(dist_innerprod(i, :));
    n1 = numel(dist_innerprod(i, :));
    mean2 = mean(dist_innerprod(A_length, :));
    var2 = var(dist_innerprod(A_length, :));  
    n2 = numel(dist_innerprod(A_length, :));
    
    metric(i) = cohen_d_unpooled(mean1, var1, n1, mean2, var2, n2);
end

% Calculate thresholds by using cutoff of D' value
thresh = xgiveny_avetemplate(cutoff, metric, A_descending);
% 10-28-23: Set default threshold to maximum stimulus level if D' is below cutoff
% cutoff at all levels. 
if isnan(thresh)
    thresh = max(A_descending);
end

% 8-27-23: 2nd criteria is that D' must monotonically increase from
% threshold to higher stimuli levels when below a certain cutoff D'. This is to avoid erroneously low
% threshold estimates from noise at lower levels (e.g. awake mouse).
willdecrease = metric(2:end) > metric(1:end-1);
willdecrease(1) = 0; % Ok for 90 dB to decrease given speaker distortion at maximum amplitude causes artifactual decrease in metric from 90 to 80 dB.
if any(willdecrease)
    ind_highest_decrease = find(willdecrease, 1);
    lowest_thresh = A_descending(ind_highest_decrease);
    
    % reset threshold to lowest level that is monotically increasing
    % even if D' value is already above cutoff.
    if metric(ind_highest_decrease) < MONOTONIC_INCREASE_BELOW_THIS_D
        thresh = max([thresh, lowest_thresh]);  % only adjust threshold if doing so raises it based on monotically increasing criteria
    end
end

%% Make plots
if ~PLOT_FIGURES
    return
end

%Plot average versus SPL
figure('DefaultAxesFontSize', 20)
plot(A_descending, metric, 'LineWidth', 3);
ylabel('Cohen D (a.u.)', 'FontSize', 24)
my_yline = yline(cutoff, '--', ['Cutoff @ D = ', num2str(cutoff)]);
title(['Threshold at ', num2str(round(thresh, 1)), ' dB'], 'FontSize', 24)
xlabel('Stimulus level (dB SPL)', 'FontSize', 24)
set(my_yline, 'FontSize', 18)
set(gca,'TickDir','out');
set(gca,'box','off')
xlim([A_descending(end)-5, A_descending(1)+5]);
ylim_bounds = [min([metric; 0]), 1.05*max(metric)];
ylim(ylim_bounds);

% Draw vertical line at threshold amplitude 
if ~isnan(thresh)
    hold on
    line([thresh thresh], [ylim_bounds(1), interp1(A_descending, metric, thresh, 'linear');], 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--'); % vertical line at exact threshold
    hold off
else
    disp('Warning: No ABR threshold detected!')
end

end
