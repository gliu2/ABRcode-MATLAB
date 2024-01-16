INDEX_INNERPROD = 10:12;
INDEX_INNERPRODWINDOW = 16:18; % rows of this metric
LEGEND_LABELS = {'8 khz pre', '16 khz pre', '32 khz pre', '8 khz post 1w', '16 khz post 1w', '32 khz post 1w'};
innerprod_window_pre = table2array(b9m6lefttable_pre(INDEX_INNERPRODWINDOW , 5));
innerprod_window_1w = table2array(b9m6lefttable_1w(INDEX_INNERPRODWINDOW , 5));
% innerprod_window_pre = table2array(b9m6lefttable_pre(INDEX_INNERPROD , 5));
% innerprod_window_1w = table2array(b9m6lefttable_1w(INDEX_INNERPROD , 5));
stimulus_levels = table2array(b9m6lefttable_pre(1, 5));

num_freq =3;
colors_p = ['r', 'g', 'b'];
figure
for ff = 1:num_freq
    plot(stimulus_levels , innerprod_window_pre(ff, :), '-o', 'color', colors_p(ff), 'LineWidth', 3)
    hold on
end
for ff = 1:num_freq
    plot(stimulus_levels , innerprod_window_1w(ff, :), '--o', 'color', colors_p(ff), 'LineWidth', 3)
    hold on
end

legend(LEGEND_LABELS, 'location', 'northwest')
set(gca,'FontSize',24)
ylabel('Inner product (nV^2)')
xlabel('Stimulus level (dB SPL)')

