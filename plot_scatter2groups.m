
function plot_scatter2groups(data, num_std, group_names, this_ylabel, dataset_names, this_ylim)

dmean = mean(data, 'omitnan');                                                 % Mean
% dci  = std(data)*tinv(0.975,size(data,1)-1);                        % Confidence Intervals
dci  = std(data, 'omitnan')*num_std;                        % Confidence Intervals
xt = [1:size(data, 2)];                                                         % X-Ticks
sep_xspace = linspace(-0.05, 0.05, 5);
xtd = repmat(xt, size(data,1), 1);                                  % X-Ticks For Data
sb = [xt'-ones(size(data,2),1)*0.1,  xt'+ones(size(data,2),1)*0.1]; % Short Bar X
lb = [xt'-ones(size(data,2),1)*0.2,  xt'+ones(size(data,2),1)*0.2]; % Long Bar X
figure('DefaultAxesFontSize', 20)
% plot(xtd + sep_xspace', data, 'o')
plot(xt, data, 'o')
hold on
for k1 = 1:size(data,2)
    plot(lb(k1,:), [1 1]*dmean(k1), sb(k1,:), [1 1]*(dmean(k1)-dci(k1)), sb(k1,:), [1 1]*(dmean(k1)+dci(k1)), '-k')
end
hold off
% set(gca, 'XTick', xt, 'XTickLabel', {'Left','Middle','Right'})
set(gca, 'XTick', xt, 'XTickLabel', group_names)
xlabel('Group')
% ylabel('Velocity (Furlongs/Fortnight)')
ylabel(this_ylabel)

legend(dataset_names, 'Location', 'bestoutside')

SPACE_SIDES = 0.5;
xlim([SPACE_SIDES, size(data, 2)+SPACE_SIDES])
ylim(this_ylim)

end