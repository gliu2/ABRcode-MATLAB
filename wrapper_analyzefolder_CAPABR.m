% function wrapper_analyzefolder_CAPABR(varargin)
% Analyze all average CAP and ABR data in folder of single trial runs of
% CAP and ABR experiments
%
% 8/3/2023 George Liu
% Dependencies: plotstack_average_CAPABR.m

%% Constants
LOAD_PATH = '\\Ricci-abr\d\George\CAP_ABR\';
SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures

% Default settings for when running function directly
FILTER_ON = 1; % Apply bandpass filter with Blackman window if no preacquisition filtering done
PLOT_FIGURES = 1; % Hard code logical for whether or not to plot figures 
SAVE_FIGURES = 1; % Hard code logical for whether or not to save figures 

selpath = uigetdir(LOAD_PATH); 
disp([' Analyzing folder ', selpath])
disp([' Filter on: ', num2str(FILTER_ON)])
disp([' Plot figures: ', num2str(PLOT_FIGURES)])
disp([' Save figures: ', num2str(SAVE_FIGURES)])

listing = dir(fullfile(selpath, '*.arf'));

allnames = {listing.name}'; % column cell array of filenames

n_filestems = length(allnames);
for i=1:n_filestems
% for i=3:n_filestems
    [~, this_filestem] = fileparts(allnames{i});
    disp(['  Analyzing filestem ', num2str(i), ' out of ', num2str(n_filestems), ': ', this_filestem])
    
    % Get first csv file with filestem
    listing2 = dir(fullfile(selpath, [this_filestem, '*.csv']));
    this_filestem_file1 = listing2(1).name;
    
    % Analyze all files with that filestem
    plotstack_average_CAPABR(selpath, this_filestem_file1, FILTER_ON, PLOT_FIGURES, SAVE_FIGURES)
end