%% analyze_innerprod_ABR.m
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
%        mylags - shifts (index, a.u.) to calculate inner products, accounting for 
%               decreasing wave latencies at higher dB levels. Time-shifts 
%               of single traces are relative to max dB average ABR trace.
%               Vector of size (A_length, 1). If matrix of size (A_length,
%               n) then returns p-values for each column vector as input.
%
%        scaleopt - normalization for xcorr (inner product) calculation,
%                   see Matlab documentation for "xcorr" for options
%                   ('none', 'biased', 'unbiased', 'coeff')
%
% Output: dist_innerprod - Inner product distributions (lag-adjusted).
%               Cell of size (A_length, 1), where A_length is # of dB levels.
%               Entry dist_innerprod{i} gives vector of size (m_traces, 1),
%               whose entry j is inner product (lag-adjusted) of j-th
%               single trace (at i-th dB level) with max dB level average ABR trace.
%
% Dependencies: none
% Last edit: 5/25/2019
%
% Author: George Liu

function dist_innerprod = analyze_innerprod_ABR(X_csv, mylags, scaleopt)
A_length = length(X_csv); % # of stimulus dB levels

% Signal basis vector: averaged ABR trace at max dB level (not normalized)
signal_basis = mean(X_csv{end}, 2); % SAMPLES x 1 vector

% Compute distribution of inner products (single traces) at each dB level
dist_innerprod = cell(A_length, 1);
for i = 1:A_length
    thisX = X_csv{i}; % SAMPLES x m_traces matrix
    samples = size(thisX, 1);
    m_traces = size(thisX, 2);
    dist_innerprod{i} = zeros(m_traces, 1); % m_traces x 1 vector
    
    % Compute each single trace inner product, using MATLAB implementation
    for j = 1:m_traces
        single_trace = thisX(:, j); % SAMPLES x 1 vector
        [r, lags] = xcorr(signal_basis, single_trace, scaleopt); % r: (2*SAMPLES-1) x 1 vector. xcorr: 2nd input shifts by lag amount relative to 1st input
        disp(size(lags))
        disp(size(mylags))
        dist_innerprod{i}(j) = r(lags'==mylags(i)); % double value
    end
end

end