%% save_allfigs.m
% 
% Save multiple open figures as PNG images
%
% Run this script after figures have been opened, to keep ordering by
% figure number
%
% Dependencies: none
% Last edit: 4/18/2019
%
% Author: George Liu

% FolderName = '/Users/georgeliu/Documents/George/KS_analyzeNoise_4-22-19';   % Your destination folder
FolderName = uigetdir;   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');

%% Save hit probs 4-18-19
for iFig = 2:length(FigList)
  FigHandle = FigList(iFig);
%   FigName   = get(FigHandle, 'Name');  % Adjust the FigName to your needs.
  FigName   = ['hitprob_', num2str(A_csv(22-iFig)), 'dB_alltimes_fp05_singletail_allnoisedist'];  % Adjust the FigName to your needs.
%   savefig(FigHandle, fullfile(FolderName, [FigName, '.png']));    %<---- 'Brackets'
  saveas(FigHandle, fullfile(FolderName, [FigName, '.png']));    % PNG file
end

%% Save hit probs 5-22-19 (from ks_2ms_noise.m cell titled, "ALTERNATIVE 4-18-19: noise dist over all traces (Find threshold versus time based on single-point analysis)")
for iFig = 2:length(FigList)
  FigHandle = FigList(iFig);
%   FigName   = get(FigHandle, 'Name');  % Adjust the FigName to your needs.
  FigName   = ['hitprob_', num2str(A_csv(22-iFig)), 'dB_alltimes_fp05_singletail_allnoisedist'];  % Adjust the FigName to your needs.
%   savefig(FigHandle, fullfile(FolderName, [FigName, '.png']));    %<---- 'Brackets'
  saveas(FigHandle, fullfile(FolderName, [FigName, '.fig']));    % FIG file
  saveas(FigHandle, fullfile(FolderName, [FigName, '.eps']));    % EPS file
  saveas(FigHandle, fullfile(FolderName, [FigName, '.png']));    % PNG file
end

%% Noise variance plots
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = ['varianceNoise_ftest_', num2str(A_csv(21-iFig)), 'dB'];  % Adjust the FigName to your needs.
  saveas(FigHandle, fullfile(FolderName, [FigName, '.png']));    % PNG file
end


%% Save 'extract_ARItissueRGB.m' script figures
FolderName = '/Users/georgeliu/Documents/Arri_analyss_4-22-19';   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');

FigNames = cell(7, 1);
FigNames{1} = [LASER, '_boxplots_4-12-19'];
FigNames{2} = [LASER, '_R_histogram_4-12-19'];
FigNames{3} = [LASER, '_G_histogram_4-12-19'];
FigNames{4} = [LASER, '_B_histogram_4-12-19'];
FigNames{5} = [LASER, '_overlayedHistogram_R_4-12-19'];
FigNames{6} = [LASER, '_overlayedHistogram_G_4-12-19'];
FigNames{7} = [LASER, '_overlayedHistogram_B_4-12-19'];

for iFig = 1:length(FigList)
  disp(['Saving figure ', num2str(iFig), ' out of ', num2str(length(FigList)), '...'])
  FigHandle = FigList(iFig);
  FigName   = FigNames{iFig};
  saveas(FigHandle, fullfile(FolderName, [FigName, '.png']));    % PNG file
end

%% 4-29-2019 Save figures using plot titles as figure file names
% FolderName = 'C:\Users\CTLab\Documents\George\Arri_analysis_4-29-19';   % Your destination folder
FolderName = 'C:\Users\CTLab\Documents\George\Arri_analysis_4-29-19_4-28data';
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 2:length(FigList)
  disp(['Saving figure ', num2str(iFig), ' out of ', num2str(length(FigList)), '...'])
  FigHandle = FigList(iFig);
  FigName   = FigHandle.Children.Title.String;
  saveas(FigHandle, fullfile(FolderName, [FigName, '.png']));    % PNG file
end