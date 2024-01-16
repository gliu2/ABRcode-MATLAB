function num = nameToNumber(x)
% Convert character array of number to int. Convert 'neg' to negative sign.
%
% Input: x - character of number, eg, '10' (-> 10), 'neg30' (-> -30)
%
% Last edit George Liu 2-11-23
% Dependencies: none

isneg = strfind(x, 'neg');
if isneg
    num = str2num(x(4:end)) * -1;
else
    num = str2num(x);
end
end