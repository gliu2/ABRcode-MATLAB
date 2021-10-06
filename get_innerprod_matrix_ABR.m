%% get_innerprod_matrix_ABR.m
% 
% Return inner product matrix of ABR single traces with all other single 
% traces, and computes average of distribution of inner products for each 
% single trace. Alternative metric to upfront inner product with coherent average. 
%
% Assumptions:
% Dataset is assumed to contain (m_traces) single-trace ABR measurements 
% with (SAMPLES) number of datapoints/timepoints
%
% Input: X - ABR single trace dataset.
%            Matrix of size (SAMPLES, m_traces), for set of single ABR
%            traces at same dB level. Columns are individual trace data.
%
% Output: innerprod_mean - vector of size (m_traces, 1),
%               whose entry j is average inner product of j-th single trace with
%               all (m_traces - 1) other single traces. 
%
% Edited 8/30/21, George Liu
%
% Dependencies: merge_singletraceABR_polarities.m
%
% Email from Daibhid 8/27/21
% Regarding an alternative to taking the inner product with the average 
% ABR. The idea is to still create a score for each ABR. Taking the inner 
% product between each single trial and all other N-1 trials will create a 
% distribution of values for that single trial. We want to generate a 
% score from that distribution. The simplest option is to take the average 
% over the N-1 inner products as the score for the single trial. Repeating 
% over all single trials will generate N scores and form the score 
% distribution. We can then analyze this distribution as usual.

function innerprod_mean  = get_innerprod_matrix_ABR(X_csv)
% % merge alternating polarities of adjacent single trace pairs to cancel out cochlear microphonic differences
% X_csv = merge_singletraceABR_polarities(X_csv);

[~, m_traces] = size(X_csv);
M = X_csv' * X_csv; % size (m_traces, m_traces)
innerprod_matrix = M - diag(diag(M));
innerprod_mean = sum(innerprod_matrix, 2)/(m_traces-1);

end

