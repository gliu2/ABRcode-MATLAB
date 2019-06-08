%% vp2p_abr_sp.m
% 
% Find peak-to-peak voltages in single traces, using same two time-points as
% were used to calculate maximum peak-to-peak voltage in average ABR trace.
%
% Assumptions:
% Peak to peak voltage is defined as difference between neighboring positive peak and
% negative peak voltage values.
%
% Input: X - ABR data (single traces).
%            Matrix of size (SAMPLES, m_traces) where columns are individual trace data.
% 
%         PEAK_THRESH (optional) - scalar for peak detection threshold. Default 0.50.
%
% Output: p2p - peak-to-peak voltages, (m_traces, 1)
%
% Dependencies: PTDetect.m
% Last edit: 6/5/2019
%
% Author: George Liu

function p2p = vp2p_abr_sp(X, varargin) 

% Identify peaks and troughs
if nargin==1
    PEAK_THRESH = 0.5;
else
    PEAK_THRESH = varargin{1};
end

% samples = size(X, 1);
m_traces = size(X, 2);

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
t1 = Pall(ind);
t2 = Tall(ind);

% Find time-aligned peak-to-peak amplitude differences in single traces
p2p = X(t1, :) - X(t2, :); % vector of length (1, m_traces)
p2p = p2p';

end