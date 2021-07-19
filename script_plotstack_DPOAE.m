% Script to plot and analyze DPOAE data 
%
% 7/12/2021 George Liu
% Dependencies: import_DPOAEcsv.m, getCDT.m, savefig_multiformat.m

%% Constants
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures

%% Load data
% [filename,path] = uigetfile('*.csv');
[M, A_csv, freq_csv, f1, f2, sample_period] = import_DPOAEcsv(filename, path);
num_samples = size(M, 2);
X = 97656.25/2048* (0:(num_samples-1));

% Analyze DPOAEs for each frequency separately
[freq_unique, ~, ic] = unique(freq_csv);
n_freq = length(freq_unique);

for ff = 1:n_freq
    disp(['Working on ', num2str(freq_unique(ff)), ' Hz...'])
    
    % amplitudes and DPOAE values for this frequency only
    A_csv2 = A_csv(ic==ff, :);
    M2 = M(ic==ff, :);
    
    num_levels = size(M2, 1);
    dB = cell(num_levels, 1);
%     dB(:) = {'dB (dBv)'};
    dB(:) = {'dB'};
    
    % Analyze DPOAE as before
    [A_descending, I] = sort(A_csv2, 'descend');
    M_sorted = M2(I, :);
    M_sorted = M_sorted / (10^9); % convert dBnv to dBv
    A_descending_cell = cellstr(num2str(A_descending));
    newYlabels=cellfun(@(x,y) [x  ' ' y], A_descending_cell, dB, 'un', 0);

    %% Plot DPOAE measurements in time domain in vertical stack
    figure
    t = tiledlayout(num_levels, 1);
    ylim_max = [-130, 10];
    for i = 1:num_levels
        nexttile
        plot(X, M_sorted(i, :), 'LineWidth', 1.5)
        ylim(ylim_max)
        ylabel(newYlabels{i})

        % Mark fundamental and cubic distortion frequencies
        [f1_out, f1_amp, f2_out, f2_amp, cdt_freq, cdt_amp, z_score] = getCDT(M_sorted(i,:));
        hold on
        scatter([f1_out, f2_out], [f1_amp, f2_amp], [], 'k', 'filled')
        if z_score > 2
            cdt_color = 'green';
        else
            cdt_color = 'red';
        end
        scatter(cdt_freq, cdt_amp, [], cdt_color, 'filled')
        hold off

        set(gca,'box','off')
        if i~=num_levels
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
        end
        
        set(gca,'FontSize',10)
    end
    t.TileSpacing = 'tight';
    t.Padding = 'tight';
    this_freq = freq_unique(ff);
    plot_title = ['DPOAE @ ', num2str(this_freq), ' Hz'];
    title(t, plot_title)
    xlabel(t, 'Frequency (Hz)')
    ylabel(t, 'Amplitude (dBv)')
    
    % Save figure
    disp('Saving figure')
    [~, save_file, ~] = fileparts(filename);
%     [save_file, save_path] = uiputfile('*.*'); % user select filepath
    save_file_freq = [save_file, num2str(this_freq), 'hz'];
    savefig_multiformat(gcf, SAVE_PATH, save_file_freq)
    
end
disp('Done')