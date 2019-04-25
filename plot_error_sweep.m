figure, plot(mrange, err_Ath_new_cache, 'LineWidth', 1, 'Color', 'red')
hold on
plot(mrange, err_Ath_old_cache, 'LineWidth', 1, 'Color', 'blue')
hold off
% xlabel('SPL_{TH}')
xlabel('Sweeps', 'FontSize', 32)
ylabel('Relative error of threshold estimate', 'FontSize', 32)
title(['Threshold estimate as function of number of sweeps (', ' {\muV}^2', ')'])
box off % remove ticks on top and right borders of plot
% axis tight % makes edges of data flush with left and right borders of plot
legend(['New method (individual traces)'],'Old method (1 averaged trace)', ...
    'location','northeast', 'FontSize', 16)
legend boxoff % remove box around legend