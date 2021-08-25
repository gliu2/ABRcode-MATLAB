function outTable = plotstack_averageABR(varargin)
% Plot a stack of average ABRs in CSV file, one stack per frequency. Saves
% plots in multiple image formats.
% Analyzes wave 1 amplitude and latency in highest stimulus waveform.
% Estimate ABR threshold.
%
% Optional input: plotstack_averageABR(path, filename)
%
% 7/11/2021 George Liu
% Dependencies: import_averageABRcsv.m, savefig_multiformat.m, get_wave1_averageABR

%% Constants
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures
SAMPLE_PERIOD = 40.96; % sample period of 40.96 us/sample in average trace

%% Load average trace data
if nargin == 2
    path = varargin{1};
    filename = varargin{2};
elseif nargin == 0
    [filename,path] = uigetfile('*.csv');
else
    disp('Warning: number of input arguments is not 0 or 2!')
    [filename,path] = uigetfile('*.csv');
end

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

% Obtain y axis scale
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

% Plot ABRs in column stacks for each frequencies, in one tiled layout
A_levels = unique(A_csv);
num_levels = length(A_levels);
fig = figure;
t = tiledlayout(num_levels, n_freq, 'TileIndexing', 'columnmajor');
for ff = 1:n_freq
    this_freq = freq_unique(ff);
%     disp(['Working on ', num2str(this_freq), ' Hz...'])
    
    % amplitudes and DPOAE values for this frequency only
    A_csv2 = A_csv(ic==ff, :);
    M2 = M(ic==ff, :);
    
%     num_levels = size(M2, 1);
    dB = cell(num_levels, 1);
    dB(:) = {'dB (nV)'};

    [A_descending, I] = sort(A_csv2, 'descend');
    M_sorted = M2(I, :);
    A_descending_cell = cellstr(num2str(A_descending));
    newYlabels=cellfun(@(x,y) [x  ' ' y], A_descending_cell, dB, 'un', 0);

    %% Plot ABRs in vertical stack
%     figure
%     t = tiledlayout(num_levels, 1);
%     s = stackedplot(X, M_sorted', 'Title', plot_title, 'DisplayLabels', newYlabels, 'LineWidth', 1.5);
    
    % Calculate threshold
    thresh_cache = get_thresh_averageABR(M_sorted, A_descending);

    % initialize variables
    amp_cache = NaN;
    lat_cache = NaN;
    is_abovenoise = zeros(num_levels, 1);
    for i = 1:num_levels
        nexttile
        y = M_sorted(i, :);
        plot(X, y, 'LineWidth', 1.5)
        ylim(ylim_max)
        
        % show ylabels for first column only
        if ff==1
            ylabel(newYlabels{i})
            % rotate y label to make horizontal
            hYLabel = get(gca,'YLabel');
            set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right')
        end

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
%         if amp > 3*noise_floor
        if A_descending(i) > thresh_cache
            text(0.99*xL(2), 0.99*yL(1), str,'HorizontalAlignment','right','VerticalAlignment','bottom', 'FontWeight', 'bold')
            is_abovenoise(i,1)=1;
        else
            text(0.99*xL(2), 0.99*yL(1), str,'HorizontalAlignment','right','VerticalAlignment','bottom')
            is_abovenoise(i,1)=0;
        end

%         % Calculate normalized cross covariance for each wave with prior level,
%         % except max level
% %         if i~=1
% %             c = xcov(y, M_sorted(i-1, :), 'normalized');
% %             num_lags = numel(c);
% %             ind_lag0 = ceil(num_lags/2);
% %             c_lag0 = c(ind_lag0);
% %             
% %             str2 = sprintf('xcor = %.2f', c_lag0);
% %             text(0.99*xL(2), 0.99*yL(2), str2,'HorizontalAlignment','right','VerticalAlignment','top')
% %         end
%         metric = get_thresh_averageABR(M_sorted);
%         str2 = sprintf('metric = %.2f', metric(i));
%         text(0.99*xL(2), 0.99*yL(2), str2,'HorizontalAlignment','right','VerticalAlignment','top')
        
        % Remove extraneous axis tick marks and x-axis from all but bottom
        % tile
        set(gca,'box','off')
        if i~=num_levels
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
        end
        
        % Remove y-axis labels from all but first column
        if ff~=1
            set(gca,'yticklabel',[])
        end
        
        % Set text size of current axis
        set(gca,'FontSize',10)
        
        % Cache wave 1 amplitude and latency for highest stimulus level
        % Show frequency as title above top-most tile of each column
        if i==1
            amp_cache = amp;
            lat_cache = lat;
            
            title([num2str(this_freq), ' Hz'])
        end
    end

    xlabel('Time (ms)')
    
%     % Save figure
% %     disp('Saving figure')
%     [~, save_file, ~] = fileparts(filename);
%     save_file_freq = [save_file, num2str(this_freq), 'hz'];
%     savefig_multiformat(gcf, SAVE_PATH, save_file_freq)
    
    % Calculate threshold based on amplitude trends
% %     ind_thresh = find(is_abovenoise, 1, 'last');
%     NOISE_MULTIPLE = 2;
%     metric = get_thresh_averageABR(M_sorted);
% %     % logistic regression to determine threshold
% %     [logitCoef,dev] = glmfit(A_descending, metric,'binomial','logit');
%     is_this_abovenoise = metric > NOISE_MULTIPLE*abs(metric(end));
%     ind_thresh = find(is_this_abovenoise, 1, 'last');
%     if isempty(ind_thresh)
%         thresh_cache = 90; % default 90 dB threshold if no waveform is above the noise floor
%     else
%         thresh_cache = str2double(A_descending_cell{ind_thresh});
%     end
    
    % Cache variables
    Amplitude(ff) = amp_cache; % nV
    Latency(ff) = lat_cache; % ms
    Threshold(ff) = thresh_cache; % dB
end

t.TileSpacing = 'none';
t.Padding = 'tight';
%     plot_title = ['ABR @ ', num2str(this_freq), ' Hz'];
%     title(t, plot_title)
%     ylabel(t, 'Amplitude (nV)')

% Maximize figure window size
fig.WindowState = 'maximized';

% Save figure
%     disp('Saving figure')
[~, save_file, ~] = fileparts(filename);
savefig_multiformat(gcf, SAVE_PATH, save_file)

% Write cached variables to excel table
[~, baseFileNameNoExt, ~] = fileparts(filename);
Filenames = cell(n_freq, 1);
Filenames(:) = {baseFileNameNoExt};  
Filenames = convertCharsToStrings(Filenames);
outTable = table(Filenames, Frequency, Amplitude, Latency, Threshold);

% % Save excel sheet output
% table_filename = [baseFileNameNoExt, '_output.xlsx'];
% writetable(outTable, fullfile(SAVE_PATH, table_filename), 'Sheet', 1, 'Range', 'A1')

% disp('Done')

end