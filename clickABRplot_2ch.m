function clickABRplot_2ch(path, exptdate, run, displayplot, saveplot, ylimitsCh1, ylimitsCh2)
% INPUTS
% path          = folder path to get to the data to plot
% exptdate      = experiment date in 'YYYY_MM_DD' format
% run           = the run number to plot
% displayplot   = toggle whether to show the plot (1 = show, 0 = don't)
% saveplot      = toggle whether or not to save a PDF (1 = save, 0 = don't)
% ylimits       = Y-axis limits (ie [-20000 5000])
%
% OUTPUTS
% Plot will be shown if displayplot = 1
% Plot will be saved in Results folder if saveplot = 1

SAVE_PATH = 'd:\users\admin\Documents\George\Results'; % path for saving figures

cd(path)
% If ylimits not commanded, use [-20 5] as default
if nargin < 6
    ylimits = [-20 5];
end

Run = run;
A = readmatrix(strcat(Run,'.csv'));

Fs = 25000;            % Sampling rate in Hz
tstep = 1/Fs;
t = [0:(length(A)-49)]; % # of time points
t = t.*tstep;           % Time vector

clf

subplot(1,2,1)
offset = 0;
for j = 2:2:19
   plot(t.*(1e3),(A(j,49:end) + offset).*0.001,'linewidth',2)
   hold on
   offset = offset - (3E3);
end
grid on
ylim(ylimitsCh1)
title(strcat(exptdate, '-', Run), 'Interpreter', 'none')
xlabel('Time (ms)')
ylabel('Voltage uV')
set(gca, 'fontsize',13)

legend('90 dB', '80 dB', '70 dB', '60 dB', '50 dB','40 dB','30 dB', '20 dB', '10 dB','location','southeast')

subplot(1,2,2)
offset = 0;
for j = 3:2:19
   plot(t.*(1e3),(A(j,49:end) + offset).*0.001,'linewidth',2)
   hold on
   offset = offset - (30E3);
end
grid on
ylim(ylimitsCh2)
title(strcat(exptdate, '-', Run), 'Interpreter', 'none')
xlabel('Time (ms)')
ylabel('Voltage uV')
set(gca, 'fontsize',13)

legend('90 dB', '80 dB', '70 dB', '60 dB', '50 dB','40 dB','30 dB', '20 dB', '10 dB','location','southeast')


if saveplot == 1
%     cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
    cd(SAVE_PATH)
    exportgraphics(gcf, strcat(exptdate,Run,'-ClickABR.pdf'), 'ContentType', 'vector');
elseif saveplot ==0
end

if displayplot == 0
    clf
elseif displayplot ==1
end

end
