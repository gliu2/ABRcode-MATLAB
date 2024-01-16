function X_csv_new = merge_singletraceABR_polarities_CM(X_csv)
% Subtract pairs of single traces in single trace ABR file
% that are alternating in polarity. Purpose is to obtain cochlear
% microphonic signal from the differences between rarefaction and condensation stimulus
% inputs.
%
% Input: X_csv - ABR single trace dataset.
%            Matrix of size (SAMPLES, m_traces), for set of single ABR
%            traces at same dB level. Columns are individual trace data.
%
% Output: X_csv_new - ABR single trace dataset, with pairs of alternating polarity subtracted and div by 2.
%            Matrix of size (SAMPLES, m_traces/2), for set of single ABR
%            traces at same dB level. Columns are individual trace data.
%
%
% Dependencies: none
% George Liu
% Last edit: 8/3/23

% If number of single traces is odd, drop last trace
if any(mod(size(X_csv, 2), 2)) % if is odd
    X_csv = X_csv(:, 1:end-1); % drop last column
end   

X_csv_new = (X_csv(:, 1:2:end) - X_csv(:, 2:2:end))/2;

end