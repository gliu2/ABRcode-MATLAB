%% extract_ARItissueRGB.m
% 
% Extract data from Ariscope RGB images of tissue, PNG or TIF.
%
% Dependencies: same_yaxes.m, same_xaxes.m
% Last edit: 4/22/2019
%
% Author: George Liu

%% Select folder containing original tissue images
% selpath = uigetdir();
% dir_csv = dir(fullfile(selpath, '*.png'));
% numfiles = length(dir_csv);

%% Select folder containing masks of tissue images
disp('Select folder containing masks of tissue images:')
selpath2 = uigetdir();
dir_csv2 = dir(fullfile(selpath2, '*.png'));
numfiles2 = length(dir_csv2);

%% Extract RGB sensor values for each tissue ROI
% if numfiles ~= numfiles2
%     disp('WARNING: Number of original images and masks not same.')
% end
% 
% filename_original = cell(numfiles, 1);
% filename_mask = cell(numfiles, 1);
% roi_r = cell(numfiles, 1);
% roi_g = cell(numfiles, 1);
% roi_b = cell(numfiles, 1);
% for i = 1:numfiles
%     disp(['Working on ', num2str(i), ' out of ', num2str(numfiles), '...'])
%     % read original tissue image PNG
%     filename_original{i} = dir_csv(i).name;
%     im1 = imread(fullfile(selpath, filename_original{i})); 
%     
%     % read corresponding mask PNG
%     filename_mask{i} = dir_csv2(i).name;
%     mask1 = imread(fullfile(selpath2, filename_mask{i})); 
%     
%     % Extract ROIs for original tissue using mask
%     mask1_bw = imbinarize(mask1);
% 
%     im1_r = im1(:,:,1);
%     im1_g = im1(:,:,2);
%     im1_b = im1(:,:,3);
% 
%     roi_r{i} = im1_r(mask1_bw(:,:,1));
%     roi_g{i} = im1_g(mask1_bw(:,:,2));
%     roi_b{i} = im1_b(mask1_bw(:,:,3));
%     
% end

%% For TIF - Extract RGB sensor values for each tissue ROI
% LASER = 'arriwhite'; % which light excitation to extract
% LASER = 'blue';
% LASER = 'green';
% LASER = 'ir7';
% LASER = 'red';
% LASER = 'violet';
% LASER = 'white17';
LASER = 'whitemix17';

disp('Select parent folder containing original tif images:')
folder_path = uigetdir();
% Get a list of all files and folders in this folder.
files = dir(folder_path);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
% Remove 'up a file' directories
% subFolders = subFolders(3:end);
subFolders(ismember( {subFolders.name}, {'.', '..'})) = [];  %remove . and ..
% % Print folder names to command window.
% for k = 1 : length(subFolders)
%   fprintf('Sub folder #%d = %s\n', k, subFolders(k).name);
% end

numfiles = length(subFolders);
if numfiles ~= numfiles2 % get numfiles2 from above for masks
    disp('WARNING: Number of original images and masks not same.')
end

filename_original = cell(numfiles, 1);
filename_mask = cell(numfiles, 1);
roi_r = cell(numfiles, 1);
roi_g = cell(numfiles, 1);
roi_b = cell(numfiles, 1);
for i = 1:numfiles
    disp(['Working on ', num2str(i), ' out of ', num2str(numfiles), '...'])
    % read original tissue image PNG
    subFolder_tif = dir(fullfile(folder_path, subFolders(i).name, 'tif', ['*', LASER,'*.tif']));
    filename_original{i} = subFolder_tif.name;
    im1 = imread(fullfile(folder_path, subFolders(i).name, 'tif', filename_original{i})); 
    
    % read corresponding mask PNG
    filename_mask{i} = dir_csv2(i).name;
    mask1 = imread(fullfile(selpath2, filename_mask{i})); 
    
    % Extract ROIs for original tissue using mask
    mask1_bw = imbinarize(mask1);

    im1_r = im1(:,:,1);
    im1_g = im1(:,:,2);
    im1_b = im1(:,:,3);

    roi_r{i} = im1_r(mask1_bw(:,:,1));
    roi_g{i} = im1_g(mask1_bw(:,:,2));
    roi_b{i} = im1_b(mask1_bw(:,:,3));
    
end

%% Plot boxplots of each tissue 'spectrum'
figure 
for j = 1:numfiles
    y2 = [roi_r{j}, roi_g{j}, roi_b{j}]; 
    g = ['R'; 'G'; 'B'];
    subplot(numfiles/3, 3, j)  
    boxplot(y2, g)
    ylabel('Intensity (a.u.)', 'FontSize', 20)
    xlabel('Sensor', 'FontSize', 20)
    title([filename_original{j}], 'Interpreter','none', 'FontSize', 20)
end

%% Plot histograms for R, G, B in 3 figures
roi_rgb = [roi_r, roi_g, roi_b];
g = ['R'; 'G'; 'B'];
for i = 1:3
    figure 
    axesHandle = zeros(numfiles, 1);
    for j = 1:numfiles
        y2 = roi_rgb{j, i};
        axesHandle(j) = subplot(numfiles, 1, j); 
        histogram(y2, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'FaceAlpha', 0.5, 'edgecolor', 'none')
        xlabel('Intensity (a.u.)')
        ylabel('Probability')
        title([filename_original{j}, ', ', g(i), ' sensor'], 'Interpreter','none')
    end
    
    % Same x and y axes for all subplots
    same_yaxes(axesHandle)
    same_xaxes(axesHandle)
end

%% Jared request: Plot histograms for R, G, B on top of each other for tissues
% create a default color map ranging from red to light pink
colors_jet = jet(numfiles); % numfiles x 3 matrix, each row is one color [r g b]

roi_rgb = [roi_r, roi_g, roi_b];
g = ['R'; 'G'; 'B'];
for i = 1:3
    figure 
    for j = 1:numfiles
        y2 = roi_rgb{j, i};
        histogram(y2, 'BinMethod', 'fd', 'Normalization', 'probability', 'LineWidth', 0.5, 'DisplayStyle', 'stairs', 'EdgeColor', colors_jet(j, :))
        if j ==1
            hold on
        end
    end
    
    xlabel('Intensity (a.u.)', 'FontSize', 28)
    ylabel('Probability', 'FontSize', 28)
    title([g(i), ' sensor'], 'Interpreter','none', 'FontSize', 28)
    
    hold off
    box off % remove ticks on top and right borders of plot
    % axis tight % makes edges of data flush with left and right borders of plot
    leg=legend(filename_original, 'location', 'northeast', 'FontSize', 14);
    set(leg,'Interpreter', 'none') % enable underscoring in legend text
    legend boxoff % remove box around legend
    
%     % Same x and y axes for all subplots
%     same_yaxes(axesHandle)
%     same_xaxes(axesHandle)
end


% %% Test individual file mask 
% [mask_file, mask_path] = uigetfile('*.png');
% mask1 = imread(fullfile(mask_path, mask_file));
% mask1_bw = im2bw(mask1);
% 
% im1_r = im1(:,:,1);
% im1_g = im1(:,:,2);
% im1_b = im1(:,:,3);
% 
% roi_r = im1_r(mask1_bw);
% roi_g = im1_g(mask1_bw);
% roi_b = im1_b(mask1_bw);
