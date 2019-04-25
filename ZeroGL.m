% Calculate zero crossing of the signal. 
% INput:  x and y. 
% Output: ZC -- position of zero crossing points
%
% Assumes unique zero crossing
% 
% 2/25/2019
% George Liu
% Dependencies: none

function ZC = ZeroGL(x,y)

y2_index = min(find(y>0)); % index of y vector just to right of 0 crossing
y1_index = max(find(y<0)); % index of y vector just to left of 0 crossing

idy = [y1_index, y2_index];

m = (y(y2_index) - y(y1_index))/(x(y2_index) - x(y1_index));
dy = 0 - y(y1_index); % from first point to zero crossing
dx = dy/m;
ZC = x(y1_index) + dx; % x-coordinate of zero crossing

% ZC = interp1(y(idy), x(idy), 0, 'linear', 'extrap' );

end

% x = linspace(0, 11*pi, 42);                                             % Create Data
% y = sin(x);                                                             % Create Data
% ZC = ZeroGL(x,y);
% figure
% plot(x, y, '-r')
% hold on
% plot(ZC, zeros(size(ZC)), 'pg')
% hold off
% grid
