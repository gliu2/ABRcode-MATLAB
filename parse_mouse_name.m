function [date, name, side] = parse_mouse_name(str)
% parse_mouse_name
% Example input: str = '20210716_b2m1_abr_left'
%
% 7/23/21 George Liu

newStr = split(str, '_');
date = newStr{1};
name = newStr{2};
side = newStr{4};

end


