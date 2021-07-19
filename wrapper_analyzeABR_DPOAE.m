% Wrapper script to run analysis on all ABR or DPOAE csv files in folder

path = uigetdir();
fileList = dir(fullfile(path, '*.csv'));
n_files = length(fileList);

for i=1:n_files
    filename = fileList(i).name;
    
    % Wrapper for analysis script
    % DPOAE
    script_plotstack_DPOAE
    
    % ABR
%     script_plotstack_averageABR
end