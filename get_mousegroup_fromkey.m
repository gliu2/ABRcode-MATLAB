function [group, date_noise] = get_mousegroup_fromkey(str_mouse_name, path_key)
% Use excel key of mouse names with noise exposure groups, to convert input
% cell array of mouse names to get groups and exposure dates.
%
% Example mouse name: b2m4 -> noise group, 20210721
%
% 7/23/21 George Liu

key = readtable(path_key);
isname = strcmp(key.Name, str_mouse_name);
group = key.Group{isname};
date_noise = key.Noise{isname};

end