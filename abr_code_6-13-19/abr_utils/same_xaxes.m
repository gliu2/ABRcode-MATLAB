%% same_xaxes.m
% 
% Set subplot x-axis scales the same
%
% Input: AxesHandles - vector of axes handles, one for each subplot
%
% Code from: https://www.mathworks.com/matlabcentral/answers/32153-subplots-with-equal-nice-y-axes-is-there-a-function
%
% Dependencies: none
% Last edit: 4/22/2019
%
% Author: George Liu

function same_xaxes(AxesHandles)

allXLim = get(AxesHandles, {'XLim'});
allXLim = cat(2, allXLim{:});
set(AxesHandles, 'XLim', [min(allXLim), max(allXLim)]);

end