%% analyze_innerprod_ABR_alllags.m
% 
% Apply inner product-sign rank method to ABR data. Returns p-value per stimulus level to
% estimate ABR threshold.
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
% Output: pval - Wilcoxin sign rank 2-sample t-test p-value, comparing
%               distribution of single-trace lag-shifted innerproducts with
%               max dB SPL average ABR trace with that of noise (0 dB)
%               single trace inner products.
%               Vector of size (A_length, 1)
%
% Dependencies: none
% Last edit: 5/25/2019
%
% Author: George Liu

function dist_innerprod = analyze_innerprod_ABR_alllags(X_csv)
A_length = length(X_csv); % # of stimulus dB levels

% Signal basis vector: averaged ABR trace at max dB level (not normalized)
signal_basis = mean(X_csv{end}, 2); % SAMPLES x 1 vector

% Compute distribution of inner products (single traces) at each dB level
dist_innerprod = cell(A_length, 1);
for i = 1:A_length
    thisX = X_csv{i}; % SAMPLES x m_traces matrix
    samples = size(thisX, 1);
    m_traces = size(thisX, 2);
    dist_innerprod{i} = zeros(m_traces, 2*samples-1); % m_traces x 2*samples-1 matrix (epoch, xcorr)
    
    % Compute each single trace inner product, using MATLAB implementation
    for j = 1:m_traces
        single_trace = thisX(:, j); % SAMPLES x 1 vector
        [r, ~] = xcorr(single_trace, signal_basis, 'coeff'); % r: (2*SAMPLES-1) x 1 vector. xcorr: 2nd input shifts by lag amount relative to 1st input
        dist_innerprod{i}(j, :) = r'; % m_traces x 2*samples-1 matrix (epoch, xcorr)
    end
end

end