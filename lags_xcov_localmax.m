%% lags_xcov_localmax.m
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

function avg_lag = lags_xcov_localmax(X_csv, dt, varargin)

% Calculate signal basis vector from normalized max averaged ABR trace
signal_basis = mean(X_csv{end}, 2); % SAMPLES x 1 vector

if nargin==2
    XCOV_THRESH_LAG = 0.5; % for coeff normalized cross-covariance
    visualizeFigs = false;
elseif nargin == 3
    XCOV_THRESH_LAG = varargin{1}; % for coeff normalized cross-covariance
    visualizeFigs = false; % turn on figures (true or false)
elseif nargin == 4
    XCOV_THRESH_LAG = varargin{1}; % for coeff normalized cross-covariance
    A_csv = varargin{2}; % turn on figures (true or false)
    visualizeFigs = false;
elseif nargin == 5
    XCOV_THRESH_LAG = varargin{1}; % for coeff normalized cross-covariance
    A_csv = varargin{2}; % turn on figures (true or false)
    visualizeFigs = varargin{3};
elseif nargin == 6
    XCOV_THRESH_LAG = varargin{1}; % for coeff normalized cross-covariance
    A_csv = varargin{2}; % turn on figures (true or false)
    visualizeFigs = varargin{3};
    signal_basis = varargin{4};
end

PEAK_THRESH = 0.50; % peak detection to find local maxima in cross-covariance plots
A_length = length(X_csv);


% Find lags using cross-correlation maxima of averaged ABR traces
avg_lag = zeros(A_length, 1);
r_cache = cell(A_length, 1);
lags_cache = cell(A_length, 1);
maxr_cache = zeros(A_length, 1);
peak_cache = cell(A_length, 1);
trough_cache = cell(A_length, 1);
for i = A_length:-1:1
    % Plot averaged ABR trace in first column
    y = mean(X_csv{i}, 2);
    [r,lags] = xcov(y, signal_basis, 'coeff');
    [P, T] = PTDetect(r, PEAK_THRESH);
    if i==A_length 
        [max_r, ind] = max(r);
        avg_lag(i) = lags(ind);
    else
        lastpeak2P = lags(P) - avg_lag(i+1);
        lag_Pnext = min(lastpeak2P(lastpeak2P>0)) + avg_lag(i+1);
        disp(i)
        
        if length(lag_Pnext)<1 % if no peak to right
            lag_Pnext = avg_lag(i+1);
        end
        
        avg_lag(i) = lag_Pnext;
        max_r = r(lags==lag_Pnext);
    end
    
    % cache vars
    r_cache{i} = r;
    lags_cache{i} = lags;
    maxr_cache(i) = max_r;
    peak_cache{i} = P;
    trough_cache{i} = T;
end

% Convert average lags to ms
avg_lag = avg_lag*dt;

if visualizeFigs
    
    % Plot xcov lags versus dB level
    hh=figure('DefaultAxesFontSize', 20);
    plot(A_csv, avg_lag, 'o-b');
    xlabel('Stimulus level (dB SPL)')
    ylabel('Lag (ms)')
    title('Average ABR cross-covariance maxima lags')
    
    % 5-25-19: Plot ave ABR cross-correlations (for finding lags) with signal basis to identify good peak
    %locations
    figure
    AxesHandles_1 = zeros(A_length, 1);
    for i = 1:A_length
        AxesHandles_1(i) = subplot(ceil(A_length/2),2,i);
        plot(lags_cache{i}*dt, r_cache{i})
        title(['Input A=', num2str(A_csv(i)), ' dB SPL'])
        xlabel('Time (ms)')
        ylabel('xcov (a.u.)')
        hold on
        % Plot peaks and troughs as circles
        scatter(lags_cache{i}(peak_cache{i})*dt, r_cache{i}(peak_cache{i}), 'b') % blue peaks
        scatter(lags_cache{i}(trough_cache{i})*dt, r_cache{i}(trough_cache{i}), 'r') % red troughs
        % Plot lag as black line
        peaklag_i = avg_lag(i);
        xline(peaklag_i, 'LineWidth', 1, 'Color', 'k'); % vertical line for cutoff
        hold off
    end
    same_yaxes(AxesHandles_1)

    % Plot maximum cross-covariance versus dB level
    figure('DefaultAxesFontSize', 20)
    plot(A_csv, maxr_cache, 'o-b')
%     hold on, yline(XCOV_THRESH_LAG); hold off % add cutoff line to max xcov vs dB plot
    xlabel('Stimulus level (dB SPL)')
    ylabel('Local max xcov (a.u.)')
    title(['Local maxima cross-covariance'])
    
end

end