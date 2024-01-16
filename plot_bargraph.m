% Bar graph
%
% Source: https://anneurai.net/2016/06/13/prettier-plots-in-matlab/
%
% Dependencies: cbrewer, ploterr, mysigstar

FONT_SIZE = 12;
THIS_XLABEL = 'Comparison group';
THIS_YLABEL = 'Difference in wave 1 amplitude measurements (nV)';

colors = cbrewer2('qual', 'Set1', 10);
% dat(:, 1) = randn(1, 50) + 10;
% dat(:, 2) = randn(1, 50) + 12;

% The data
v = wave1amp_avg;
nv = numel(v); % Number of elements
dv = abs(bsxfun(@minus,v,v')); % Absolute pairwise differences
pair_dif_wave1avg = triu(dv, 1); % Obtain vector of differences
pair_dif_wave1avg = pair_dif_wave1avg(pair_dif_wave1avg~=0);

dat1 = pair_dif_wave1avg;
dat2 = wave1amp_dif;
% remove nan values
dat1 = dat1(~isnan(dat1));
dat2 = dat2(~isnan(dat2));
dat = {dat1, dat2};

CATEGORIES = {['Inter-subject (n=', num2str(length(dat1)), ')'], ['Intra-subject (n=', num2str(length(dat2)), ')']};
N_CATEGORIES = length(CATEGORIES);

% Plot
figure;
% subplot(4,7,6); % rather than a square plot, make it thinner
hold on;
% if we want each bar to have a different color, loop
for b = 1:N_CATEGORIES
    bar(b, mean(dat{b}), 'FaceColor',  colors(b, : ), 'EdgeColor', 'none', 'BarWidth', 0.6);
end
 
% show standard deviation on top
h = ploterr(1:2, [mean(dat1), mean(dat2)], [], [std(dat1), std(dat2)], 'k.', 'abshhxy', 0);
set(h(1), 'marker', 'none'); % remove marker
 
% label what we're seeing
% if labels are too long to fit, use the xticklabelrotation with about -30
% to rotate them so they're readable
set(gca, 'xtick', [1 2], 'xticklabel', CATEGORIES, ...
    'xlim', [0.5 2.5]);
ylabel(THIS_YLABEL); 
xlabel(THIS_XLABEL);
 
% if these data are paired, show the differences
% plot(dat', '.k-', 'linewidth', 0.2, 'markersize', 2);
 
% significance star for the difference
% [~, pval] = ttest(dat1, dat2);
n1 = length(dat1);
n2 = sum(~isnan(dat2));
[tprob_2tail, tprob_1tail] = ttest2_mean_ste(n1, mean(dat1), std(dat1)/sqrt(n1), n2, mean(dat2), std(dat2)/sqrt(n2));
pval = tprob_2tail;
% if mysigstar gets 2 xpos inputs, it will draw a line between them and the
% sigstars on top
mysigstar(gca, [1 2], 170, pval);
 
% add significance stars for each bar
for b = 1:2
    [~, pval] = ttest(dat{b});
    yval = mean(dat{b}) * 0.5; % plot this on top of the bar
    mysigstar(gca, b, yval, pval);
    % if mysigstar gets just 1 xpos input, it will only plot stars
end

set(gca,'FontSize', FONT_SIZE)