%% - 2023_5_20
clf
clear
exptdate = '2023_5_19';
% folder = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/";
folder = "//Ricci-abr/d/Julien/";
path = strcat(folder,exptdate,'/');

%% - Click ABRs
runs = 1; %[1 29 52 7 8 24 25 28 32 33 34 49 50];
for j = 1:length(runs)
%run = 'Run1';
run = strcat('Run',num2str(runs(j)));
displayplot = 1;
saveplot = 1;
clickABRplot_2ch(path, exptdate, run, displayplot, saveplot, [-25 6],[-250 60])
end