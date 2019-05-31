%% fourier_innerprod_ABR.m
% 
% Return inner product distributions of single traces with average trace 
% ``template" at max dB in ABR data. 
%
% Assumptions:
% Dataset is assumed to contain (m_traces) single-trace ABR measurements 
% with (SAMPLES) number of datapoints/timepoints, measured at (A_length) 
% stimulus dB SPL levels. 
%
% Input: X - ABR dataset.
%            Cell of size (A_length, 1), where A_length is # of dB levels. 
%            Each cell entry X{i} gives matrix of size (SAMPLES, m_traces), 
%            where columns are individual trace data.
%
%        Fs - sampling rate (Hz)
%
%        signal_basis (optional) - signal basis ``template" in time domain
%
% Output: f - frequency domain (Hz) for plotting on x-axis. Row vector.
%
%         dist_innerprod - Inner product distributions computed in Fourier domain.
%               Cell of size (A_length, 1), where A_length is # of dB levels.
%               Entry dist_innerprod{i} gives vector of size (m_traces, 1),
%               whose entry j is inner product (lag-adjusted) of j-th
%               single trace (at i-th dB level) with max dB level average ABR trace.
%
% Dependencies: fft_p1.m
% Last edit: 5/30/2019
%
% Author: George Liu

function [f, dist_innerprod] = fourier_innerprod_ABR(X_csv, Fs, varargin)

% Signal basis vector: averaged ABR trace at max dB level (not normalized)
if nargin == 3
    signal_basis = varargin{1}; % e.g. for chunk signal bases
else
    signal_basis = mean(X_csv{end}, 2); % SAMPLES x 1 vector
end
[f, P1_basis] = fft_p1(signal_basis, Fs); % f:row, P1_basis:column vectors of length (SAMPLES/2+1)

% Compute distribution of inner products (single traces) at each dB level
A_length = length(X_csv); % # of stimulus dB levels
dist_innerprod = cell(A_length, 1);
for i = 1:A_length
    thisX = X_csv{i}; % SAMPLES x m_traces matrix
    [~, this_P1] = fft_p1(thisX, Fs); % matrix of size (SAMPLES/2+1) x m_traces

    % Compute inner product (coeff normalization) in Fourier domain of single-sided spectra
    denom = vecnorm(this_P1).*vecnorm(P1_basis); % row vector of length m_traces
    dist_innerprod{i} = this_P1' * P1_basis ./ denom'; % m_traces x 1 vector
%     dist_innerprod{i} = this_P1' * P1_basis;
end

end