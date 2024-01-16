function A = load_humanABR_singletracedata(varargin)
%Script for loading and visualizing single trial human ABR data
%
% Optional inputs:
%   1. Filename (full path). 
%   2. High pass filter (Hz)
%   3. Low pass filter (Hz)
%
% Output: 
%   A - Single trace data is loaded as a M x N matrix. Each column is a single
%       trace vector. M=1024 (samples per trace), N = 500 (# single traces).
%
% Last edit George Liu 2-15-23
% Dependencies: none

FILE_NAME = 'D:\users\admin\Documents\George\Human Single Trial ABR\SINGS001_tiph_80blk.TXT';
SAMPLE_RATE = 40000; % Hz; sampling rate is 40000 Hz, or sample period of 25 us

if nargin == 0
    filename = FILE_NAME;
elseif nargin == 1
    filename = varargin{1};
elseif nargin == 2
    filename = varargin{1};
    highpass_freq = varargin{2};
elseif nargin == 3
    filename = varargin{1};
    highpass_freq = varargin{2};
    lowpass_freq = varargin{3};
elseif nargin == 4
    filename = varargin{1};
    highpass_freq = varargin{2};
    lowpass_freq = varargin{3};
    notch_freq = varargin{4};
end

A = readmatrix(filename);

if nargin == 2
    A = highpass(A, highpass_freq, SAMPLE_RATE);
elseif nargin == 3
    A = highpass(A, highpass_freq, SAMPLE_RATE);
    A = lowpass(A, lowpass_freq, SAMPLE_RATE);
elseif nargin == 4
    A = highpass(A, highpass_freq, SAMPLE_RATE);
    A = lowpass(A, lowpass_freq, SAMPLE_RATE);
    d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1', notch_freq - 1, 'HalfPowerFrequency2', notch_freq + 1, ...
               'DesignMethod','butter','SampleRate', SAMPLE_RATE);
    A = filtfilt(d, A);       
end

end