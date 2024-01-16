% Script to plot power spectrogram of average trace
filename = 'D:\users\admin\Documents\George\Human Single Trial ABR\SINGS001_tiph_neg30blk.TXT';
Fs = 40000;

single_traces = load_humanABR_singletracedata(filename);
filtered_single_traces = load_humanABR_singletracedata(filename, 600, 3000, 60);

openLoop = mean(single_traces, 2);
buttLoop = mean(filtered_single_traces, 2);

[popen,fopen] = periodogram(openLoop,[],[],Fs);
[pbutt,fbutt] = periodogram(buttLoop,[],[],Fs);

figure
plot(fopen,20*log10(abs(popen)),fbutt,20*log10(abs(pbutt)),'--')
ylabel('Power/frequency (dB/Hz)')
xlabel('Frequency (Hz)')
title('Power Spectrum')
legend('Unfiltered','Filtered')
grid

xlim([0, 3100])