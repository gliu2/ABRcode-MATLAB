%% jasaFigures_summary_6_12_19.m
% 
% Create figures for JASA paper on inner product calculation of ABR
% single-trace ensemble threshold
%
% Run this script after "import_ABRcsv_folder.m", which loads single-trace 
% ABR data.
%
% Dependencies: same_yaxes.m, same_xaxes.m, PTDetect.m, 
%               analyze_innerprod_ABR.m, innerprod2pval.m, 
%               lags_xcov.m, vp2p_abr.m, vp2p_abr_sp.m, xgiveny.m,
%               vp2p_abr_sp_loc.m
%               
% Last edit: 6/12/2019
%
% Author: George Liu

% load carsmall MPG              % the sample dataset variable
% MPG(:,2)=MPG(:,1).*2;
% MPG(:,3)=MPG(:,1).*3;

NUM_FEATURES = 2;
NUM_DATASETS = 5;

% Select path to Excel data
disp('Select Excel file containing summary stats... ')
[file,path] = uigetfile('*.xlsx');
disp(fullfile(path, file))

%% Load excel data
dataCorners = "D3:M18";
M = xlsread(fullfile(path, file), dataCorners);

% get feature labels for each column
[~,feature_names,~] = xlsread(fullfile(path, file), 'D2:M2');

% get dataset name labels for each column
[~,dataset_names,~] = xlsread(fullfile(path, file), 'D1:M1');
for i = 2:2:length(dataset_names)+1
dataset_names{i} = dataset_names{i-1};
end
dataset_names_sparse = cell(length(dataset_names)/2, 1);
for i = 1:length(dataset_names_sparse)
    dataset_names_sparse{i} = dataset_names{1+2*(i-1)};
end

% get stat y-axis labels for each row
[~,stat_ylabel,~] = xlsread(fullfile(path, file), 'B3:B18');

% get stat summary names for each row
[~,stat_names,~] = xlsread(fullfile(path, file), 'C3:C18');

% dataset colors
dataset_colors = jet(NUM_DATASETS);

%% Collect data into groups
ylim_all = {...
    [0, 0.06], ...
    [0, 0.5], ...
    [0, 95], ...
    [0, 95], ...
    [0, 95], ...
    [0, 95], ...
    [0, 95], ...
    [0, 1], ...
    [0, 1], ...
    [0, 1], ...
    [0, 1], ...
    [0.5, 1], ...
    [0.5, 1], ...
    [0.5, 1], ...
    [0.5, 1], ...
    [0.5, 1], ...
    };

num_summarystats = length(M);
vpp_stats = M(:, 1:2:end);
ip_stats = M(:, 2:2:end);

for i=1:num_summarystats
    this_vppstat = vpp_stats(i, :);
    this_ipstat = ip_stats(i, :);
    
    data = [this_vppstat', this_ipstat'];
    plot_scatter2groups(data, 1, feature_names(1:2), stat_ylabel(i), dataset_names_sparse, ylim_all{i});
    title(stat_names{i})
    
%     figure('DefaultAxesFontSize', 20)
%     boxplot(data, 'Notch','on','Labels', feature_names(1:2), 'Whisker',1)
%     title(stat_names{i})
%     ylabel(stat_ylabel{i})
%     
% %     lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
% %     set(lines, 'Color', 'g');
% %     % Change the boxplot color from blue to green
% %     a = get(get(gca,'children'),'children');   % Get the handles of all the objects
% %     %t = get(a,'tag');   % List the names of all the objects 
% %     %box1 = a(7);   % The 7th object is the first box
% %     set(a, 'Color', 'r');   % Set the color of the first box to green
%     hold on
%     x0=ones(length(data)).*(1+(rand(length(data))-0.5)/5);
%     x1=ones(length(data)).*(1+(rand(length(data))-0.5)/10);
%     
%     MARKER_SIZE = 36; % default
%     f1=scatter(x0(:,1), data(:,1), MARKER_SIZE, dataset_colors, 'filled');f1.MarkerFaceAlpha = 0.4;hold on 
%     f2=scatter(x1(:,2).*2, data(:,2), MARKER_SIZE, dataset_colors,'filled');f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on
end

%%
disp('Select folder for saving figures... ')
FolderName = uigetdir;   % Your destination folder
disp(FolderName)

FigList = findobj(allchild(0), 'flat', 'Type', 'figure');

%Plot
DATASET_DATE_IDENTIFIER = '_6-12-19';

FIGNAMES = stat_names;

% for iFig = 1:length(FigList)
for i = 1:length(FigList)
    iFig = length(FigList) - i + 1;
    disp(['Saving figure ', num2str(i), ' out of ', num2str(length(FigList))])
    FigHandle = FigList(iFig);
    FigName   = FIGNAMES{i};  % Adjust the FigName to your needs.
    saveas(FigHandle, fullfile(FolderName, [FigName, DATASET_DATE_IDENTIFIER, '.fig']));    % FIG file
    saveas(FigHandle, fullfile(FolderName, [FigName, DATASET_DATE_IDENTIFIER, '.eps']));    % EPS file
end

disp('Done.')