%% lags_xcov.m
% 
% Find latency lags of averaged ABR traces using cross-covariance with coeff normalization
%
% Assumptions:
% Dataset is assumed to contain (m_traces) single-trace ABR measurements 
% with (SAMPLES) number of datapoints/timepoints, measured at (A_length) 
% stimulus dB SPL levels. 
%
% Input: X_csv - ABR dataset.
%            Cell of size (A_length, 1), where A_length is # of dB levels. 
%            Each cell entry X{i} gives matrix of size (SAMPLES, m_traces), 
%            where columns are individual trace data.
% 
%         XCOV_THRESH_LAG (optional) - scalar for minimum cross-covariance
%            maximum for beginning linear extrapolation of lags. Default 0.50.
% 
%         visualizeFigs (optional) - boolean for showing figures
%
% Output: avg_lag_xcovExtrap - lags, (A_length, 1) vector of lag in ms for
%             each dB level
%
% Dependencies: same_yaxes.m
% Last edit: 5/29/2019
%
% Author: George Liu

function avg_lag_xcovExtrap = lags_xcov(X_csv, dt, varargin)

if nargin==2
    XCOV_THRESH_LAG = 0.5; % for coeff normalized cross-covariance
    visualizeFigs = false;
elseif nargin == 3
    XCOV_THRESH_LAG = varargin{1}; % for coeff normalized cross-covariance
    visualizeFigs = false; % turn on figures (true or false)
elseif nargin == 4
    XCOV_THRESH_LAG = varargin{1}; % for coeff normalized cross-covariance
    A_csv = varargin{2}; % turn on figures (true or false)
    visualizeFigs = true;
end

A_length = length(X_csv);

% Calculate signal basis vector from normalized max averaged ABR trace
signal_basis = mean(X_csv{end}, 2); % SAMPLES x 1 vector

% Find lags using cross-correlation maxima of averaged ABR traces
avg_lag = zeros(A_length, 1);
r_cache = cell(A_length, 1);
lags_cache = cell(A_length, 1);
maxr_cache = zeros(A_length, 1);
for i = 1:A_length
    % Plot averaged ABR trace in first column
    y = mean(X_csv{i}, 2);
    [r,lags] = xcov(y, signal_basis, 'coeff');
    [max_r, ind] = max(r);
    avg_lag(i) = lags(ind)*dt;
    
    % cache vars
    r_cache{i} = r;
    lags_cache{i} = lags;
    maxr_cache(i) = max_r;
end

% Extrapolate starting from lowest dB where max cross-covariance is > 0.50
% Find lowest dB
[row, ~] = find(maxr_cache > XCOV_THRESH_LAG);
ind_uselag = min(row);

% Extrapolate
avg_lag_extrap = avg_lag(ind_uselag:end);
avg_lag_extrap = interp1(avg_lag_extrap, -(ind_uselag-2):1, 'linear', 'extrap')';

%TODO 5-27-19: Ensure that extrapolated lags have ceiling of max possible lag

% Save final extrapolation-correct xcov lags
avg_lag_xcovExtrap = avg_lag;
avg_lag_xcovExtrap(1:ind_uselag-1) = avg_lag_extrap(1:ind_uselag-1);


if visualizeFigs
    
    % Plot xcov lags versus dB level
    hh=figure('DefaultAxesFontSize', 20);
    plot(A_csv, avg_lag, 'o-b');
    xlabel('Stimulus level (dB SPL)')
    ylabel('Lag (ms)')
    title('Average ABR cross-covariance maxima lags')

    % Overlay plot of extrapolation on lags vs dB
    figure(hh)
    hold on
    plot(A_csv(1:length(avg_lag_extrap)), avg_lag_extrap, 'o--g')
    title(['Cross-covariance lags, max(xcov)>', num2str(XCOV_THRESH_LAG)])
    hold off
    
    % 5-25-19: Plot ave ABR cross-correlations (for finding lags) with signal basis to identify good peak
    %locations
    figure
    AxesHandles_1 = zeros(A_length, 1);
    for i = 1:A_length
        AxesHandles_1(i) = subplot(ceil(A_length/2),2,i);
        plot(lags_cache{i}, r_cache{i})
        title(['Input A=', num2str(A_csv(i)), ' dB SPL'])
        xlabel('Time (ms)')
        ylabel('xcov (a.u.)')
    end
    same_yaxes(AxesHandles_1)

    % Plot maximum cross-covariance versus dB level
    figure('DefaultAxesFontSize', 20)
    plot(A_csv, maxr_cache, 'o-b')
    hold on, yline(XCOV_THRESH_LAG); hold off % add cutoff line to max xcov vs dB plot
    xlabel('Stimulus level (dB SPL)')
    ylabel('Max xcov (a.u.)')
    title(['Max cross-covariance, cutoff=', num2str(XCOV_THRESH_LAG)])
    
end

end