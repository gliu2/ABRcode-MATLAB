function visualize_innerprod_hist(X_csv, x, A_csv, dist_innerprod)

A_length = size(X_csv, 1);

% Plot averaged ABR trace and distribution of single-trace signal
% components per dB SPL level
count2 = 0;
figure('DefaultAxesFontSize', 16)
AxesHandles_3 = zeros(A_length, 1);
axesHandle = zeros(A_length, 1);
for i = 1:A_length
    % Plot averaged ABR trace in first column - absolute scale
    count2 = count2 + 1;
    AxesHandles_3(i) = subplot(A_length, 2, count2);  
    y = mean(X_csv{i}, 2);
    plot(x, y*10^6)
    title(['Average ABR, Input A=', num2str(A_csv(i)), ' dB SPL'])
    xlabel('Time (ms)')
    ylabel('\muV')
    
    % Plot histogram of single-trace inner products
    count2 = count2 + 1;
    axesHandle(i) = subplot(A_length, 2, count2); 
    histogram(dist_innerprod{i}, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'edgecolor', 'none')
    xlabel('Inner product (a.u.)')
    ylabel('Prob')
    title(['Input A=', num2str(A_csv(i)), ' dB SPL'])

    % plot mean as blue line
    hold on
    signal_mean = mean(dist_innerprod{i});
    line([signal_mean signal_mean], ylim, 'LineWidth', 1, 'Color', 'b'); % vertical line for cutoff
    hold off
end
% Same y axes for average ABR plots
same_yaxes(AxesHandles_3)
% Same x and y axes for all histograms
same_yaxes(axesHandle)
same_xaxes(axesHandle)

end