function output_table = plotstack_averageABR(filename, path)
% Plot a stack of average ABRs in CSV file 
% Analyze wave 1 amplitude and latency in highest stimulus waveform.
% Also estimate ABR threshold.
%
% 7/11/2021 George Liu
% Dependencies: import_averageABRcsv.m, savefig_multiformat.m, get_wave1_averageABR

%% Constants
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures
SAMPLE_PERIOD = 40.96; % sample period of 40.96 us/sample in average trace

%% Load average trace data
% [filename,path] = uigetfile('*.csv');
disp(['Opening ', fullfile(path, filename)])
[M, A_csv, freq_csv] = import_averageABRcsv(filename, path);
num_samples = size(M, 2);
X = SAMPLE_PERIOD / 1000 * (1:num_samples); 

% Analyze DPOAEs for each frequency separately
[freq_unique, ~, ic] = unique(freq_csv);
n_freq = length(freq_unique);

% Initialize cache variables
Frequency = freq_unique;
Amplitude = zeros(n_freq, 1);
Latency = zeros(n_freq, 1);
Threshold = zeros(n_freq, 1);

for ff = 1:n_freq
    this_freq = freq_unique(ff);
    disp(['Working on ', num2str(this_freq), ' Hz...'])
    
    % amplitudes and DPOAE values for this frequency only
    A_csv2 = A_csv(ic==ff, :);
    M2 = M(ic==ff, :);
    
    num_levels = size(M2, 1);
    dB = cell(num_levels, 1);
    dB(:) = {'dB (nV)'};

    [A_descending, I] = sort(A_csv2, 'descend');
    M_sorted = M2(I, :);
    A_descending_cell = cellstr(num2str(A_descending));
    newYlabels=cellfun(@(x,y) [x  ' ' y], A_descending_cell, dB, 'un', 0);

    %% Plot ABRs in vertical stack
    figure
    t = tiledlayout(num_levels, 1);
%     s = stackedplot(X, M_sorted', 'Title', plot_title, 'DisplayLabels', newYlabels, 'LineWidth', 1.5);

    ylim_max = [-1200, 1200];
    % Check to make sure bounds of plot are within y limits
    ymax_data = max(M(:));
    ymin_data = min(M(:));
    if ymax_data > ylim_max(2)
        ylim_max(2) = ymax_data;
    end
    if ymin_data < ylim_max(1)
        ylim_max(1) = ymin_data;
    end
    
    % initialize variables
    amp_cache = NaN;
    lat_cache = NaN;
    thresh_cache = NaN;
    is_abovenoise = zeros(num_levels, 1);
    for i = 1:num_levels
        nexttile
        y = M_sorted(i, :);
        plot(X, y, 'LineWidth', 1.5)
        ylim(ylim_max)
        ylabel(newYlabels{i})
        
        % rotate y label to make horizontal
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right')

        % Mark wave 1 peak and following trough
        [peak_pt, trough_pt, amp, lat] = get_wave1_averageABR(X, y);
        hold on
        scatter(peak_pt(1), peak_pt(2), [], 'green') % peak
        scatter(trough_pt(1), trough_pt(2), [], 'red') % trough
        hold off
        
        % Display wave 1 measurements in corner of plot
        xL=xlim;
        yL=ylim;
        str = sprintf('P-P_I = %.0f nV', amp);
        noise_floor = std(M_sorted(end, :));
        if amp > 3*noise_floor
            text(0.99*xL(2), 0.99*yL(1), str,'HorizontalAlignment','right','VerticalAlignment','bottom', 'FontWeight', 'bold')
            is_abovenoise(i,1)=1;
        else
            text(0.99*xL(2), 0.99*yL(1), str,'HorizontalAlignment','right','VerticalAlignment','bottom')
            is_abovenoise(i,1)=0;
        end

%         % Calculate cross correlation for each wave with prior level,
%         % except max level
%         if i~=1
%             [c, ~] = xcorr(y, M_sorted(i-1, :), 50, 'normalized');
%             c_max = max(c);
%             
%             str2 = sprintf('xcor = %.2f', c_max);
%             text(0.99*xL(2), 0.99*yL(2), str2,'HorizontalAlignment','right','VerticalAlignment','top')
%         end
        
        % Remove extraneous axis tick marks and x-axis from all but bottom
        % tile
        set(gca,'box','off')
        if i~=num_levels
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
        end
        
        set(gca,'FontSize',10)
        
        if i==1
            amp_cache = amp;
            lat_cache = lat;
        end
    end
    t.TileSpacing = 'tight';
    t.Padding = 'tight';
    plot_title = ['ABR @ ', num2str(this_freq), ' Hz'];
    title(t, plot_title)
%     ylabel(t, 'Amplitude (nV)')
    xlabel('Time (ms)')
    
    % Save figure
    disp('Saving figure')
    [~, save_file, ~] = fileparts(filename);
    save_file_freq = [save_file, num2str(this_freq), 'hz'];
    savefig_multiformat(gcf, SAVE_PATH, save_file_freq)
    
    % Calculate threshold based on amplitude trends
    ind_thresh = find(is_abovenoise, 1, 'last');
    thresh_cache = str2double(A_descending_cell{ind_thresh});
    % Cache variables
    Amplitude(ff) = amp_cache; % nV
    Latency(ff) = lat_cache; % ms
    Threshold(ff) = thresh_cache; % dB
end

% Write cached variables to excel table
outTable = table(Frequency, Amplitude, Latency, Threshold);
[~, baseFileNameNoExt, ~] = fileparts(filename);
table_filename = [baseFileNameNoExt, '_output.xlsx'];
writetable(outTable, fullfile(SAVE_PATH, table_filename), 'Sheet', 1, 'Range', 'A1')

disp('Done')

end