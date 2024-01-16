function [peak_pt, trough_pt, amp, latency, latency_trough] = get_wave1_humanaverageABR(x, y, varargin)
% Measure wave 1 in average ABR waveform using peak detection algorithm.
% Calculates peak-peak wave 1 amplitude from wave 1 peak to following
% trough.
%
% Intput: 
%       x - vector of x-coordinates of average ABR plot
%       y - vector of y-coordinates of average ABR plot
%       prior_wave1_latency (optional)        - prior wave 1 peak latency at higher stimulus level
%       prior_wave1_trough_latency (optional) - prior wave 1 trough latency at higher stimulus level
%
% Output: 
%       peak_pt - 2 element vector of wave 1 peak coordinates, (x1,y1)
%       trough_pt - 2 element vector of wave 1 trough coordinates, (x2,y2)
%       amp     - scalar value of wave 1 p-p amplitude. Equal to y1-y2.
%       latency - scalar value of wave 1 latency. Equal to x1.
%       latency_trough
%
% Last edit George Liu
% 2/12/23 - adapted from get_wave1_averageABR.m. Lower PEAK_SENSITIVITY
%   from 0.1 to 0.01 to capture more subtle troughs
% 12/8/21 - if trough is not found, use prior trough latency to estimate
% trough location
%
% 10/22/21 - added latency of trough as an output
%
% 8/31/2021 
%       Updated to add optional input to specify x-coordinate of wave 1
%       peak, and rules to ensure wave 1 peak and trough latencies do not
%       decrease at lower stimulus levels.
%
% Dependencies: PTDetect.m

%% Constants
PEAK_SENSITIVITY = 0.03;
% Parameters for searching for wave 1 peak in highest stimulus level
PEAK_START = 1; % start time (ms) 
PEAK_END = 3.5;
TROUGH_START = 1.5; % ms
TROUGH_END = 4; % consider changing to 2.9
% Parameters for searching for wave 1 peak in highest stimulus level
MAX_PEAK_TROUGH_TIME = 1; % ms. Maximum time between wave 1 peak and trough. Consider changing to 0.9
MAX_PEAK_LAT_BACKWARD = 0; % ms. Maximum time that wave 1 peak latency can decrease at lower stimulus level.

%% get wave 1 peak coordinates
[P, T] = PTDetect(y, PEAK_SENSITIVITY);

x_peak = x(P);
y_peak = y(P);

is_ind_peak_candidates = (x_peak > PEAK_START) & (x_peak < PEAK_END);
if any(is_ind_peak_candidates)
    [peak_y, ~] = max(y_peak(is_ind_peak_candidates));
    ind = find(y_peak == peak_y);
    peak_pt = [x_peak(ind), y_peak(ind)];
else
    is_ind_peak_candidates_alt = (x > PEAK_START) & (x < PEAK_END);
    [peak_y, ~] = max(y(is_ind_peak_candidates_alt));
    ind = find(y == peak_y);
    peak_pt = [x(ind), y(ind)];
end

% get wave 1 following trough coordinates
x_trough = x(T);
y_trough = y(T);

is_ind_trough_candidates = (x_trough > TROUGH_START) & (x_trough < TROUGH_END) & (x_trough > peak_pt(1));
if ~any(is_ind_trough_candidates)
    if isempty(varargin)
        disp('Warning: Trough not found. In get_wave1_humanaverageABR (line 67')
        peak_pt = [NaN, NaN];
        trough_pt = [NaN, NaN];
        amp = peak_pt(2) - trough_pt(2);
        latency = peak_pt(1);
        latency_trough = trough_pt(1);
        return
    else
        % 12-8-21: if trough not found, then estimate trough using prior wave 1 latency
        disp('Warning: Trough not found; using previous level trough to estimate. In get_wave1_averageABR (line 75')
        prior_wave1_trough_latency = varargin{2};
        
        % find candidate trough point closest to prior trough
        [~, ind_x_closest_trough] = min(abs(x_trough - prior_wave1_trough_latency));
        trough_pt = [x_trough(ind_x_closest_trough), y_trough(ind_x_closest_trough)];
    end
else
    ind2 = find(is_ind_trough_candidates, 2);
    trough_pt = [x_trough(ind2(1)), y_trough(ind2(1))];
    % Check that next trough is not real trough
    if length(ind2) > 1
        if (y_trough(ind2(2)) < y_trough(ind2(1))) && (x_trough(ind2(2)) < TROUGH_END)
            trough_pt = [x_trough(ind2(2)), y_trough(ind2(2))];
        end
    end
end

if ~isempty(varargin)
    prior_wave1_latency = varargin{1};
    while (trough_pt(1) < peak_pt(1)) || (peak_pt(1) < prior_wave1_latency - MAX_PEAK_LAT_BACKWARD)
        
    % Added 8-31-21
    % If input includes optional wave 1 peak latency of higher stimulus level,
    % then ensure wave 1 peak latency does not decrease
        
        % Edit wave 1 latency only if it has decreased from higher stimulus
        % level
        if peak_pt(1) < prior_wave1_latency - MAX_PEAK_LAT_BACKWARD
            % determine trough first, as wave 1 peak must be before
            trough_max_time = prior_wave1_latency + MAX_PEAK_TROUGH_TIME;
            is_ind_trough_candidates_next = (x_trough > prior_wave1_latency) & (x_trough < trough_max_time);

            % If trough exists in time window after prior wave 1
            % peak latency, use it.
            if any(is_ind_trough_candidates_next)
                % identify trough with smallest value
                y_trough_candidate = inf;
                ind_trough_candidate = find(is_ind_trough_candidates_next);
                n_trough_candidates = numel(ind_trough_candidate);
                for i=1:n_trough_candidates
                    this_y_trough_candidate = y_trough(ind_trough_candidate(i));
                    if this_y_trough_candidate < y_trough_candidate
                        y_trough_candidate = this_y_trough_candidate;
                        x_trough_candidate = x_trough(ind_trough_candidate(i));
                    end
                end
            else
            % If no trough in time window after prior wave 1 peak latency, then
            % estimate trough latency using prior latency (if optional input is
            % provided). Otherwise use local minimum in time window.
                if numel(varargin)==2
                    % Use prior wave 1 latency if not before 
                    prior_wave1_trough_latency = varargin{2};
                    ind_trough_candidate = find(x==prior_wave1_trough_latency, 2);
                    y_trough_candidate = y(ind_trough_candidate(1));
                    x_trough_candidate = x(ind_trough_candidate(1));

                    % Check that next trough is not real trough
                    if length(ind_trough_candidate) > 1
                        if y(ind_trough_candidate(2)) < y(ind_trough_candidate(1))
                            y_trough_candidate = y(ind_trough_candidate(2));
                            x_trough_candidate = x(ind_trough_candidate(2));
                        end
                    end
                else
                    % use local minimum in time window
                    ind_prior_wave1_latency = find(x==prior_wave1_latency, 1);
                    ind_trough_max_time = find(x==trough_max_time, 1);
                    trough_search_window = ind_prior_wave1_latency:ind_trough_max_time;
                    [y_trough_candidate, ind3] = min(y(trough_search_window));
                    x_trough_candidate = x(trough_search_window(ind3));
                end
            end

            % Update new trough point
            trough_pt = [x_trough_candidate, y_trough_candidate];

            % Determine new wave 1 peak, now that trough is determined. 
            is_ind_peak_candidates_next = (x_peak > prior_wave1_latency) & (x_peak < trough_pt(1));
            if ~any(is_ind_peak_candidates_next)
                % Keep same wave 1 latency if no peaks identified
                peak_pt(1) = prior_wave1_latency;
                peak_pt(2) = y(x==prior_wave1_latency);
            else
                % Find local peak with largest amplitude
                ind_peak_candidates_next = find(is_ind_peak_candidates_next);
                n_peak_candidates = numel(ind_peak_candidates_next);
                new_peak_candidate = -inf;
                for i=1:n_peak_candidates
                    this_new_peak_candidate = y_peak(ind_peak_candidates_next(i));
                    if this_new_peak_candidate > new_peak_candidate
                        new_peak_candidate = this_new_peak_candidate;
                        x_peak_candidate = x_peak(ind_peak_candidates_next(i));
                    end
                end
                peak_pt(1) = x_peak_candidate;
                peak_pt(2) = new_peak_candidate;
            end
        end
    

        % Added 2/12/23
        % Ensure that wave 1 trough follows (and does not come before) wave 1 peak
        if trough_pt(1) < peak_pt(1)
            % Use first candidate trough after peak
            is_trough_after_peak = x_trough > peak_pt(1);
            ind_troughs_after_peak = find(is_trough_after_peak);
            if any(is_trough_after_peak)
                first_trough_after_peak = x_trough(ind_troughs_after_peak(1));
                trough_pt(1) = first_trough_after_peak;
                trough_pt(2) = y_trough(ind_troughs_after_peak(1));
            end
        end
    end
end

% calculate wave 1 amplitude and latency
amp = peak_pt(2) - trough_pt(2);
latency = peak_pt(1);
latency_trough = trough_pt(1);

end