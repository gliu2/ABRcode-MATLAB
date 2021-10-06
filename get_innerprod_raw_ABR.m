%% get_innerprod_raw_ABR.m
% 
% Return inner product distribution of a single traces set with their average trace 
% ``template" at the same dB in ABR data. 
%
% Inner products are calculated as unbiased cross-variance (lag = 0) between a single trace ABR
% at a given dB level and the coherent average ABR at same dB level. 
%
% Cross-covariance is calculated as raw (see revision note at bottom of comments,except changed to raw). 
%
% Assumptions:
% Dataset is assumed to contain (m_traces) single-trace ABR measurements 
% with (SAMPLES) number of datapoints/timepoints
%
% Input: X - ABR single trace dataset.
%            Matrix of size (SAMPLES, m_traces), for set of single ABR
%            traces at same dB level. Columns are individual trace data.
%
% Output: innerprod - vector of size (m_traces, 1),
%               whose entry j is inner product of j-th single trace with
%               same dB level coherent average ABR trace.
%
% Note: adapted from analyze_innerprod_ABR.m, which used coherent average
% from maximum dB level instead of from same dB level.
%
% Dependencies: merge_singletraceABR_polarities.m
% Last edit: 8/30/2021
%
% Author: George Liu
%
% Included revision to analyze_innerprod_ABR_DOM.m by Daibhid O Maoileidigh, Last edit: 7/23/2019
% Change the template to be the average at each SPL level
% To account for larger values at larger SPLs, use the unbiased cross correlation

function innerprod  = get_innerprod_raw_ABR(X_csv)
% % merge alternating polarities of adjacent single trace pairs to cancel out cochlear microphonic differences
% X_csv = merge_singletraceABR_polarities(X_csv); 

% [numpts, m_traces] = size(X_csv);
X_avg = mean(X_csv, 2); % size (numpts, 1) vector
% X_csv_norm = diag(sqrt(X_csv' * X_csv)); % size (m_traces, 1) vector
% X_avg_norm = sqrt(X_avg' * X_avg); % double

% save inner product with various normalization schemes
innerprod = X_csv' * X_avg; % size (m_traces, 1) vector. Equivalent to raw cross correlation with no lag. 
% innerprod_norm = X_csv' * (X_avg./X_avg_norm); % size (m_traces, 1) vector. 
% innerprod_norm2 = (X_csv'./X_csv_norm) * (X_avg./X_avg_norm); % size (m_traces, 1) vector. Equivalent to normalized cross correlation.
end
