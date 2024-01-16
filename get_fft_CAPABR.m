function [f, P1, phi] = get_fft_CAPABR(y_avg)
% Calculate single-Sided Amplitude Spectrum of single trace (or average
% trace) CAP/ABR FFT.
%
% Assumes sampling rate used in acquiring single trial trace.
%
% Can plot output using:
%     plot(f,P1,'LineWidth', 3)
%
% 7/12/2023 George Liu
% Dependencies: none

SAMPLE_PERIOD = 5.12; % sample period of 5.12 us/sample in single trace
SAMPLE_PERIOD_MS = SAMPLE_PERIOD/1000; % sample period in ms
Fs = 1000/SAMPLE_PERIOD_MS; % Sampling frequency = 195312.50 Hz

y_fft = fft(y_avg);
L = length(y_avg); % sample length 
P2 = abs(y_fft/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

phi = angle(y_fft/L); % units of [-pi, pi]
phi = phi(1:L/2+1);

% plot(f,P1,'LineWidth', 3) 
% % FFT plot labels
% title("Single-Sided Amplitude Spectrum of X(t)")
% xlabel("f (Hz)")
% ylabel("|P1(f)|")

end