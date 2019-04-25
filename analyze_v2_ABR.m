%% analayze_v2_ABR.m
% 
% Analyze ABR signal for a single tone burst
%
% Input: X - ABR signals, SAMPLES x m matrix - where m ABR experiments each of SAMPLES length
% Output: V2 - mean square voltage across time of ABR response(s), m x 1
%               vector
%
% Dependencies: none
% Last edit: 2/6/2019
%
% Author: George Liu

function V2 = analyze_v2_ABR(X)
V2 = diag(X'*X)/size(X,1);

end
