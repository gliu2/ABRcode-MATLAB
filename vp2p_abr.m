%% vp2p_abr.m
% 
% Find maximum peak-to-peak voltage(s) in ABR trace(s)
%
% Assumptions:
% Peak to peak voltage is defined as difference between neighboring positive peak and
% negative peak voltage values.
%
% Input: X - ABR data.
%            Matrix of size (SAMPLES, m_traces) where columns are individual trace data.
% 
%         PEAK_THRESH (optional) - scalar for peak detection threshold. Default 0.50.
%
% Output: max_p2p - maximum peak-to-peak voltages, (m_traces, 1)
%
% Dependencies: PTDetect.m
% Last edit: 7/27/23
%
% Author: George Liu

function max_p2p = vp2p_abr(X, varargin) 

% Identify peaks and troughs
if nargin==1
    PEAK_THRESH = 0.5;
else
    PEAK_THRESH = varargin{1};
end

% samples = size(X, 1);
m_traces = size(X, 2);

max_p2p = zeros(m_traces, 1);
for i = 1:m_traces
    y = X(:, i); % single trace ABR; volts (V)
    y_rel = y/max(abs(y)); % y vector is scaled to max 1

    [P, T] = PTDetect(y_rel, PEAK_THRESH);

    %Find maximum peak-to-peak amplitude in average ABR at each dB level
    % If no peak, then throw warning and return NAN
    if isempty(P) || isempty(T)
        max_p2p(i) = NaN;
        disp(['  Warning: no peak found in trace ', num2str(i)])
        continue
    end
    
    % If peak first, then pair peak with trough and trough before
    if P(1) < T(1) % # troughs = # peaks OR # peaks-1
        peak2peak_pt = y(P(1:length(T))) - y(T); % peak then trough
        peak2peak_tp = y(P(2:end)) - y(T(1:length(P)-1)); % trough then peak
    else % T(1) < P(1) % # troughs = # peaks OR # peaks+1
        peak2peak_pt = y(P(1:length(T)-1)) - y(T(2:end));
        peak2peak_tp = y(P) - y(T(1:length(P)));
    end
    
    % Find maximum peak-to-peak amplitude
    [max_p2p(i), ~] = max([peak2peak_pt; peak2peak_tp]);

end

end