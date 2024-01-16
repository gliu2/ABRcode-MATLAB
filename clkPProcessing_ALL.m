%Script for grand average click + data
%
% Last edit George Liu 1-22-23
% Dependencies: none

clear all
close all
%%

% directory1 = 'C:\Users\gmusacch\Documents\MATLAB\EPProcessing\FFR\';
% directory1 = 'G:\ALL_Data\SEI\clkplus\';
directory1 = 'D:\users\admin\Documents\George\Human Single Trial ABR';

f = dir(fullfile(directory1, '*.txt'));

% f = dir([ directory1 '*100*.txt']);
% f2 = dir([ directory1 '*200*.txt']);
% f3 = dir([ directory1 '*ah*.txt']);
% f4 = dir([ directory1 '*clk*80*.txt']);


%FILENAME GENERATION
for i=1:1:size(f,1)
    ff = f(i).name;
    filenames{i} = ff;
end
% for i=1:1:size(f2,1)
%     ff2 = f2(i).name;
%     filenames2{i} = ff2;
% end
% % 
% for i=1:1:size(f3,1)
%     ff3 = f3(i).name;
%     filenames3{i} = ff3;
% end
% for i=1:1:size(f4,1)
%     ff4 = f4(i).name;
%     filenames4{i} = ff4;
% end


%% MAKE da GNDs
% j=3; %columns of data from excel file, 3 corresponds to average amplitude value converted to uV
% k=2;
% d1 = designfilt('lowpassiir','FilterOrder',24, ...
%     'HalfPowerFrequency',0.9,'DesignMethod','butter');

j=2; %column where average uV data reside 
k=2; %row where first data point resides
time=-15.9:0.0250:60.88; 

for ii=1:size(filenames,2)
    temp = dlmread(fullfile(directory1, filenames{ii}), ',', 21,1);
    m=0; %%number of cols between avg uV data conditions
    for iii=1:(size(temp,2)-1)/6
        resp = temp(1:3072,j+m);
%         time = temp(1:3072,2);
        %     respflt=filtfilt(d1,resp);
        %     clear b,     clear a,
        adrate = 40000;
        ckpALL(ii,:,iii) = resp;
        clear resp
        m=m+6;
    end
end


%% PLOT individual time waveforms
ckpSUB=squeeze(ckpALL(3,:,:));

h=figure;
subplot(5,1,1)
plot(time, ckpSUB(:,[1,2,7,11])); box 'off'
xlim([-5 60]); ylim([-1 1])
xlabel ('Time(ms)'); ylabel('Amplitude')
legend ('clk1', 'clk2', 'clk3', 'clk4');legend('boxoff');
title('CKP005')

subplot(5,1,2)
plot(time, ckpSUB(:,[3,4])); box 'off'
xlim([-5 60]); ylim([-1 1])
xlabel ('Time(ms)'); ylabel('Amplitude')
legend ('clkAH100', 'AH100 alone');legend('boxoff');

subplot(5,1,3)
plot(time, ckpSUB(:,[11,10])); box 'off'
xlim([-5 60]); ylim([-1 1])
xlabel ('Time(ms)'); ylabel('Amplitude')
legend ('clkAH200', 'AH200 alone');legend('boxoff');

subplot(5,1,4)
plot(time, ckpSUB(:,[5,6])); box 'off'
xlim([-5 60]); ylim([-1 1])
xlabel ('Time(ms)'); ylabel('Amplitude')
legend ('clkSAW100', 'SAW100 alone');legend('boxoff');

subplot(5,1,5)
plot(time, ckpSUB(:,[9,8])); box 'off'
xlim([-5 60]); ylim([-1 1])
xlabel ('Time(ms)'); ylabel('Amplitude')
legend ('clkSAW200', 'SAW200 alone');legend('boxoff');


%% FFT
Fs = 40104;           % Sampling frequency of IHS &USB lite (3072 pts/76.78 ms)*1000 =40104 pts/sec
T = 1/Fs;             % Sampling period

L = 1391;   % Length of signal (after clipping to FFR) CAPD
t = (0:L-1)*T;        % Time vector
npow = 6;
n = 2^nextpow2(L)*npow;

%For clk+
f0beg=1437; %20 ms
f0end=3072; %end

DA=ckpSUB;
% DA=clkOneHun;
% DA=clkAh;
for iv = 1:size(ckpSUB,2)
    X = ckpSUB(f0beg:f0end,iv);
    Y = fft(X, n);
    P2 = abs(Y/n);
    P1 = P2(1:n/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(n/2))/n;
    ckpSUBfft(:,iv) = P1*npow; %need to multiply the amplitude by the factor that is padded in order to preserve correct amps
    %     clear X, clear Y
    clear P1, clear P2
end



%% PLOT
f_fft=f;

h=figure;
subplot(5,1,1)
plot(f_fft, ckpSUBfft(:,[1,2,7,11])); box 'off'
xlim([0 1000]); ylim([0 .2])
xlabel('f (Hz)'), ylabel('|P1(f)|')
legend ('clk1', 'clk2', 'clk3', 'clk4');legend('boxoff');
title('FFT CKP005')

subplot(5,1,2)
plot(f_fft, ckpSUBfft(:,[3,4])); box 'off'
xlim([0 1000]); ylim([0 .2])
xlabel('f (Hz)'), ylabel('|P1(f)|')
legend ('clkAH100', 'AH100 alone');legend('boxoff');

subplot(5,1,3)
plot(f_fft, ckpSUBfft(:,[11,10])); box 'off'
xlim([0 1000]); ylim([0 .2])
xlabel('f (Hz)'), ylabel('|P1(f)|')
legend ('clkAH200', 'AH200 alone');legend('boxoff');

subplot(5,1,4)
plot(f_fft, ckpSUBfft(:,[5,6])); box 'off'
xlim([0 1000]); ylim([0 .2])
xlabel('f (Hz)'), ylabel('|P1(f)|')
legend ('clkSAW100', 'SAW100 alone');legend('boxoff');

subplot(5,1,5)
plot(f_fft, ckpSUBfft(:,[9,8])); box 'off'
xlim([0 1000]); ylim([0 .2])
xlabel('f (Hz)'), ylabel('|P1(f)|')
legend ('clkSAW200', 'SAW200 alone');legend('boxoff');


%% F0, HFs
% TAKE TOTAL RMS OF THE FFR PERIOD

% f0range = (24:48); % F0, 75-175 bins 9-19 Skoe 2015
% f1range = (19:78); % F1 - 175-750  19-78
% allrange = (9:78);



f0range = (24:55); % 100 Hs
f1range = (55:72); % 200 Hz
allrange = (72:247); % HFs

% DAfft = clkAhfft;
DAfft = clkTwoHunfft;
DAfft = clkOneHunfft;
for vi = 1:size(DAfft,2)
    DAfftAmps(vi,1) = rms(DAfft(f0range, vi)); %RMS of F0
    DAfftAmps(vi,2) = rms(DAfft(f1range, vi)); %RMS of F1
    DAfftAmps(vi,3) = rms(DAfft(allrange, vi)); %RMS of whole range
end
clear vi


%%
figure
subplot(3,1,1)
bar(clkAhfftAmps); box('off')
ylim([0 .09])
legend ('75-150 Hz', '175-225 Hz', '225-800 Hz'); legend('boxoff');
ylabel('|P1(f)|')
title('Mean RMS of Ah stimuli')

subplot(3,1,2)
bar(clkOneHunfftAmps);box('off')
legend ('75-150 Hz', '175-225 Hz', '225-800 Hz'); legend('boxoff');
ylim([0 .09])
ylabel('|P1(f)|')
title('Mean RMS of 100 Hz stimuli')


subplot(3,1,3)
bar(clkTwoHunfftAmps); box('off')
legend ('75-150 Hz', '175-225 Hz', '225-800 Hz'); legend('boxoff');
ylim([0 .09])
ylabel('|P1(f)|')
title('Mean RMS of 200 Hz stimuli')


