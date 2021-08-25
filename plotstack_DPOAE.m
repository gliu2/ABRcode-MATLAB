function outTable = plotstack_DPOAE(varargin)
% Plot and analyze DPOAE data in one CSV file.
%
% 7/12/2021 George Liu
% Dependencies: import_DPOAEcsv.m, getCDT.m, savefig_multiformat.m

%% Constants
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures

%% Load data
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

[M, A_csv, freq_csv, f1, f2, ~] = import_DPOAEcsv(filename, path);
num_samples = size(M, 2);
X = 97656.25/2048* (0:(num_samples-1));

% Analyze DPOAEs for each frequency separately
[freq_unique, ~, ic] = unique(freq_csv);
n_freq = length(freq_unique);
num_levels = length(unique(A_csv));

% Initialize cache variables
Frequency = freq_unique;
Threshold = zeros(n_freq, 1);

% Plot each frequency in a separate column in one tiled layout
fig = figure;
t = tiledlayout(num_levels, n_freq, 'TileIndexing', 'columnmajor');
ylim_max = [-130, 10];
for ff = 1:n_freq
%     disp(['Working on ', num2str(freq_unique(ff)), ' Hz...'])
    
    % amplitudes and DPOAE values for this frequency only
    A_csv2 = A_csv(ic==ff, :);
    M2 = M(ic==ff, :);
    this_f1 = f1(ic==ff);
    this_f2 = f2(ic==ff);
    
    dB = cell(num_levels, 1);
%     dB(:) = {'dB (dBv)'};
    dB(:) = {'dB (dBv)'};
    
    % Analyze DPOAE as before
    [A_descending, I] = sort(A_csv2, 'descend');
    M_sorted = M2(I, :);
    M_sorted = M_sorted / (10^9); % convert dBnv to dBv
    A_descending_cell = cellstr(num2str(A_descending));
    newYlabels=cellfun(@(x,y) [x  ' ' y], A_descending_cell, dB, 'un', 0);

    %% Plot DPOAE measurements in time domain in vertical stack
    is_abovenoise = zeros(num_levels, 1);
    for i = 1:num_levels
        nexttile
        plot(X, M_sorted(i, :), 'LineWidth', 1.5)
        ylim(ylim_max)
        
        % show y label of first column only
        if ff==1
            ylabel(newYlabels{i})
            
            % rotate y label to make horizontal
            hYLabel = get(gca,'YLabel');
            set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right')
        end

        % Mark fundamental and cubic distortion frequencies
        [f1_out, f1_amp, f2_out, f2_amp, cdt_freq, cdt_amp, z_score] = getCDT(M_sorted(i,:));
%         [f1_out, f1_amp, f2_out, f2_amp, cdt_freq, cdt_amp, z_score] = getCDT(M_sorted(i,:), this_f1(i), this_f2(i));
        hold on
        scatter([f1_out, f2_out], [f1_amp, f2_amp], [], 'k', 'filled')
        if z_score > 2
            is_abovenoise(i) = 1;
            cdt_color = 'green';
        else
            is_abovenoise(i) = 0;
            cdt_color = 'red';
        end
        scatter(cdt_freq, cdt_amp, [], cdt_color, 'filled')
        hold off
        
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
        
        set(gca,'FontSize',10)
        
        % title of top tile of each column
        if i==1
            this_freq = freq_unique(ff);
            title([num2str(this_freq), ' Hz'])
        end
        
        % show xlabel at bottom of each column
        if i==num_levels
            xlabel('Frequency (Hz)')
        end
    end
    
    % Calculate threshold based on amplitude trends
    ind_thresh = find(is_abovenoise, 1, 'last');
    if isempty(ind_thresh)
        thresh_cache = 90; % default 90 dB threshold if no waveform is above the noise floor
    else
        thresh_cache = str2double(A_descending_cell{ind_thresh});
    end
    % Cache variables
    Threshold(ff) = thresh_cache; % dB
end

t.TileSpacing = 'none';
t.Padding = 'tight';
% this_freq = freq_unique(ff);
% plot_title = ['DPOAE @ ', num2str(this_freq), ' Hz'];
% title(t, plot_title)
% xlabel(t, 'Frequency (Hz)')
% ylabel(t, 'Amplitude (dBv)')

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
outTable = table(Filenames, Frequency, Threshold);

% disp('Done')

end