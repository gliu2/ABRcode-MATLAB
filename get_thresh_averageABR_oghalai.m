function [thresh, metric] = get_thresh_averageABR_oghalai(M_sorted, A_descending, varargin)
% Measures maximum peak-peak amplitude of average ABR trace as metric to 
% determine ABR treshold, when metric is 4 standard deviations above the
% noise floor, which is the RMS of the average noise trace. 
%
% Input: M_sorted - sorted average ABR traces. Matrix of size (n_traces x n_samples), 
%                   each row is average ABR trace at level, in
%                   descending intensity level. All same frequency.
%        A_descending - vector of stimulus levels (dB) in descending order
%        PLOT_FIGURES (optional) - boolean, specify whether to plot figures
%
% output: 
%       thresh - threshold level (dB) calculated using normalized cross
%                covariance between average ABR traces ('coeff') at zero lag. 
%       metric - normalized cross covariance values between average ABR
%                traces.
%
% Method as originally described by Oghalai for offline CAP threshold analysis:
%   John S. Oghalai, "Chlorpromazine inhibits cochlear function in guinea
%   pigs," Hearing Research, 198:1–2 (59-68), 2004.
%   https://doi.org/10.1016/j.heares.2004.03.013.
%   https://www.sciencedirect.com/science/article/pii/S0378595504001212
%
% Later Oghalai/Ricci papers that use for offline ABR threshold analysis:
% (1) Wenzel et al. "Laser-induced collagen remodeling and deposition 
% within the basilar membrane of the mouse cochlea". J Biomed Opt (2007).
% (2) Xia et al. "Deficient forward transduction and enhanced reverse transduction
% in the alpha tectorin C1509G human hearing loss mutation", Disease Models
% and Mechanisms (2010). 
% (3)* Cho et al. "Mechanisms of Hearing Loss after Blast Injury to the
% Ear". PLOS ONE (2013).
% (4)* Becker et al. "The presynaptic ribbon maintains vesicle populations 
% at the hair cell afferent fiber synapse". eLIFE (2018).
%
% Note: adapted code from jasaFigures_6_8_19_DOM.m
%
% Last edit: 8/31/21 George Liu
%
% Dependencies: vp2p_abr.m, xgiveny_avetemplate.m

%% Constant
PEAK_THRESH = 0.1;
PLOT_FIGURES = true;

if nargin==3
    PLOT_FIGURES = varargin{1};
end

%% Calculate ABR threshold and level function using max peak-peak amplitude
%Find maximum peak-to-peak amplitude in average ABR at each dB level
A_length = length(A_descending);
max_p2p = zeros(A_length, 1);
for i=1:A_length
    max_p2p(i) = vp2p_abr(M_sorted(i,:)', PEAK_THRESH);
end

% Calculate threshold based on noise floor
noise_rms = rms(M_sorted(end, :)); % RMS of noise floor

% Calculate standard deviation of noise signal (lowest stimulus level). 
% Could approximate with RMS unless signal average is not zero.
% noise_std = std(M_sorted(end, :));
noise_std_alllevels = std(M_sorted'); % row vector containing standard deviation for each level
noise_std = min(noise_std_alllevels); % in case there is extra noise in 0 dB trace c/w 10 dB trace

% Calculate thresholds by fitting ABR p-p voltages with a cubic spline and 
% identifying when signal is four standard deviations above the noise floor
p2p_cutoff = noise_rms + 4*noise_std;
A_query = A_descending(end):A_descending(1);
p2p_spline_fit = interp1(A_descending, max_p2p, A_query, 'spline');
thresh = xgiveny_avetemplate(p2p_cutoff, p2p_spline_fit, A_query);
metric = max_p2p;

%% Make plots
if ~PLOT_FIGURES
    return
end

%Plot average versus SPL
figure('DefaultAxesFontSize', 20)
plot(A_query, p2p_spline_fit, 'b');
ylabel('Max peak-to-peak amplitude (nV)', 'FontSize', 32)
my_yline = yline(p2p_cutoff, '--', ['THRESHOLD = ', num2str(p2p_cutoff), ' nV']);
title(['Peak-to-peak amplitude of average ABR'], 'FontSize', 32)
xlabel('Stimulus level (dB SPL)', 'FontSize', 32)
set(my_yline, 'FontSize', 18)
set(gca,'TickDir','out');
set(gca,'box','off')
xlim([A_descending(end)-5, A_descending(1)+5]);
ylim_bounds = [min([p2p_spline_fit, 0]), 1.05*max(p2p_spline_fit)];
ylim(ylim_bounds);

% Draw vertical line at threshold amplitude 
if ~isnan(thresh)
    hold on
    line([thresh thresh], [ylim_bounds(1), interp1(A_descending, max_p2p, thresh, 'spline');], 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--'); % vertical line at exact threshold
    hold off
else
    disp('Warning: No ABR threshold detected by Oghalai peak-peak method!')
end

end
