%% find_xgiveny.m
% 
% Find x-coordinate where plot first reaches input y-value.
% Uses linear interpolation.
% Returns NaN if no intercept (x-value for input y-value)
%
% For ABR threshold given cutoff (feature y-axis, dB x-axis)
%
% Input:  cutoff - y-value cutoff for finding 1st x-intercept
%
%         y_features - Feature data per dB level.
%            Vector of size (A_length, 1)
%
%         A_csv - dB levels / x-axis values for plotting y.
%
% Output: x_thresh - x-coordinate where plot first reaches input y-value
%
% Dependencies: none
% Last edit: 6/7/2019
%
% Author: George Liu

function x_thresh = xgiveny(y_cutoff, y_features, A_csv)

is_thresh = y_features > y_cutoff;
if any(is_thresh) && ~all(is_thresh) % ensure curve is not entirely below or above cutoff
    Athresh_ind = find(is_thresh, 1) - 1;
    
    if Athresh_ind == 0
        x_thresh = NaN;
    else
%         Athresh_ind = max(Athresh_ind, 1); % Ensure if curve above cutoff, then threshold is lowest value
        A_thresh = A_csv(Athresh_ind);

        % Linear interpolate to identify exact amplitude corresponding to
        % 3*RMS_0
        dify = y_features(Athresh_ind + 1) - y_features(Athresh_ind);
        difx = A_csv(Athresh_ind + 1) - A_thresh; % should be 5 dB or whatever amplitude spacing is
        mslope = dify/difx;
        deltay = y_cutoff - y_features(Athresh_ind);
        x_thresh = A_thresh + deltay/mslope;
    end
    
else
%     disp('Warning: No xgiveny threshold!')
    x_thresh = NaN;
end

end