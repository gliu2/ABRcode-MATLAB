function [peak_pt, trough_pt, amp, lat_peak, lat_trough] = get_wave1_averageCAP(X, y_avg, varargin)
% Calculate wave 1 amplitude for CAP. Flip CAP to make sure wave 1 peak is positive.
% CAP traces are inverted compared to ABR. Then analyze as would for ABR.
%
% Last edit: 8/26/23 George Liu
%
% Dependencies: get_wave1_averageABR.m

if nargin==2
    [peak_pt, trough_pt, amp, lat_peak, lat_trough] = get_wave1_averageABR(X, -1*y_avg);
elseif nargin==3
    [peak_pt, trough_pt, amp, lat_peak, lat_trough] = get_wave1_averageABR(X, -1*y_avg, varargin{1});
elseif nargin==4
    [peak_pt, trough_pt, amp, lat_peak, lat_trough] = get_wave1_averageABR(X, -1*y_avg, varargin{1}, varargin{2});
else
    error('Invalid number of input arguments to get_wave1_averageCAP')
end
peak_pt(2) = -1 * peak_pt(2);
trough_pt(2) = -1 * trough_pt(2);

end