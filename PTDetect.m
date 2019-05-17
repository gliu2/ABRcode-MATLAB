function [P,T] = PTDetect(x,E)
% Usage [P,T] = PTDetect(x,E), in which E is a threshold, x is a vector of data, P is the list of peaks, and T is the list of troughs.
% Based on Jacobson 2001 "Auto-threshold Peak Detection in Physiological
% Signals ", Papers from 23rd Annual International Conference of the IEEE Engineering in Medicine and Biology Society, October
% 25-28, 2001, held in Istanbul, Turkey. See also ADM001351 for entire conference on cd-rom., The original document
% contains color images.
P = [];
T = [];
a = 1;
b = 1;
i = 0;
d = 0;
xL = length(x);
 while i~= xL
    i = i+1;
    if d==0
        if x(a) >= x(i)+E
            d = 2;
        elseif x(i) >= x(b)+E
            d = 1;
        end
        if x(a) <= x(i)
            a = i;
        elseif x(i) <= x(b)
            b = i;
        end
    elseif d==1
        if x(a) <= x(i)
            a = i;
        elseif x(a) >= x(i)+E
            P = [P a]; b = i; d = 2;
        end
    elseif d==2
        if x(i) <= x(b)
            b = i;
        elseif x(i) >= x(b)+E
            T = [T b]; a = i; d = 1;
        end
    end
 end
end