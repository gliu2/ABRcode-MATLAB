function [date, name, studytype, side, metadata] = get_mousefile_metadata(filename)
% Get metadata of ABR file based on filename. 
% Includes date of ABR experiment, mouse label, and side of ABR.
%
% Input: filename - name of file. String. Optional to include extension or
% full path.
%
% Output:
%       date - date of study, e.g., '20210829'. String / char array
%
% Last edit 12/7/21, George Liu

%% Constants
PATH_KEY = 'd:\users\admin\Documents\George\mouse_metadata.xlsx';

%% Breakdown file name  
[~, baseFileNameNoExt, ~] = fileparts(filename);
filename_words = split(baseFileNameNoExt, "_");
date = convertCharsToStrings(filename_words{1});
name = filename_words{2};
studytype = convertCharsToStrings(filename_words{3});
side = convertCharsToStrings(filename_words{4});

%% Obtain metadata from key file
metadata_table = readtable(PATH_KEY);

is_mouse = strcmpi(name,metadata_table.Mouse_name);
metadata = metadata_table(is_mouse, :);
% disp(metadata)

% delta_time = daysact('7-sep-2002',  '25-dec-2002')

end