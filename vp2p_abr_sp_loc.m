%% vp2p_abr_sp_loc.m
% 
% Find locations of the maximum peak and minimum trough in coherent average
% ABR trace that were used to calculate maximum peak-to-peak voltage.
%
% Assumptions:
% Peak to peak voltage is defined as difference between neighboring positive peak and
% negative peak voltage values.
%
% Input: X - ABR data (single traces) at one dB level.
%            Matrix of size (SAMPLES, m_traces) where columns are individual trace data.
% 
%         PEAK_THRESH (optional) - scalar for peak detection threshold. Default 0.50.
%
% Output: ind_peak - index of peak for peak-peak voltage in coherent average BR
%         ind_trough - index of trough for peak-peak voltage in coherent average BR
%
% Dependencies: PTDetect.m
% Last edit: 6/8/2019
%
% Author: George Liu

function [ind_peak, ind_trough] = vp2p_abr_sp_loc(X, varargin) 

% Identify peaks and troughs
if nargin==1
    PEAK_THRESH = 0.5;
else
    PEAK_THRESH = varargin{1};
end

% Find two time points of peak and trough in average ABR that give its max
% peak-peak voltage
y = mean(X, 2); % SAMPLES x 1 vector
y_rel = y/max(abs(y)); % y vector is scaled to max 1

[P, T] = PTDetect(y_rel, PEAK_THRESH);

%Find maximum peak-to-peak amplitude in average ABR

% If peak first, then pair peak with trough and trough before
if P(1) < T(1) % # troughs = # peaks OR # peaks-1
    peak2peak_pt = y(P(1:length(T))) - y(T); % peak then trough
    peak2peak_tp = y(P(2:end)) - y(T(1:length(P)-1)); % trough then peak
    Pall = [P(1:length(T)), P(2:end)];
    Tall = [T, T(1:length(P)-1)];
else % T(1) < P(1) % # troughs = # peaks OR # peaks+1
    peak2peak_pt = y(P(1:length(T)-1)) - y(T(2:end));
    peak2peak_tp = y(P) - y(T(1:length(P)));
    Pall = [P(1:length(T)-1), P];
    Tall = [T(2:end), T(1:length(P))];
end

% Find maximum peak-to-peak amplitude
[~, ind] = max([peak2peak_pt; peak2peak_tp]);
ind_peak = Pall(ind);
ind_trough = Tall(ind);

end