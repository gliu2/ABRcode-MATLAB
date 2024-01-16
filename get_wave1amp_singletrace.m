function dist = get_wave1amp_singletrace(X_csv, timewindow_peak_range, timewindow_trough_range)
% Calculate distribution of wave 1 amplitude values for single trace ABRs,
% using average of multiple points to calculate peak and trough values in
% single traces prior to taking their difference.
%
% 12-29-21

X_csv_peak = X_csv(timewindow_peak_range, :);
X_csv_trough = X_csv(timewindow_trough_range, :);

dist = mean(X_csv_peak, 1) - mean(X_csv_trough, 1);
dist = reshape(dist, [length(dist), 1]); % make column vector