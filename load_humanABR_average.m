function [x, y, y_uv] = load_humanABR_average(varargin)
%Script for loading and visualizing average trace human ABR data
%
% Output:
%   x - vector of time, units of ms
%   y - vector of amplitudes, units of uV
%
% Last edit George Liu 2-11-23
% Dependencies: none
FILE_NAME = 'D:\users\admin\Documents\George\Human Single Trial ABR\SINGS001_tiph_80avg.TXT';
% TIME_UNIT = 'ms';
% VOLTAGE_UNIT = 'uV';


if nargin == 0
    filename = FILE_NAME;
else
    filename = varargin{1};
end

A = readmatrix(filename); % Automatically starts importing data from index 1 of average human ABR file

% Check that first row is index 1
if A(1,1) ~= 1
    error(['Automatic parameter detection mis-loaded average human ABR data for file ', filename])
end

x = A(:, 2); % time (ms)
y = A(:, 3); % voltage (A.U.)
y_uv = A(:, 4); % voltage (uV)

end