% Merge single trace human ABRs of opposite polarities (C - condensation
% and R - rarefaction).
%
% Hard coded 'path' variable must contain single trace ABR data files
% (.txt.) with the following file name format:
%   ST02_tipv_-10_C_st.TXT
%
% 1-1-2024
%
% Dependencies: load_humanABR_singletracedata.m

% SINGLE TRACE KEYWORD
SINGLETRACE_KEY = 'st'; % single trace (which can be condensation, rarefaction, or bulk average). Unlike bulk average which is already the average of rarefaction and condensation traces

%% Hard code these parameters
% MONTAGE = 'tiph';
MONTAGE = 'tipv';
% MONTAGE = 'wick';
% path = 'D:\George-abr\human_abr\SingleTrials_aas_01\Matt Initial Recordings';
% path = 'D:\George-abr\human_abr\SingleTrials_aas_01\Matt Repeat Recordings_st';
% path = 'D:\George-abr\human_abr\SingleTrials_aas_02\Izzy Initial Recordings_st';
% path = 'D:\George-abr\human_abr\SingleTrials_aas_02\Izzy Repeat Recordings_st';
% path = 'D:\George-abr\human_abr\SingleTrials_aas_03\George Initial Recordings_st';
% path = 'D:\George-abr\human_abr\SingleTrials_aas_03\George Repeat Recordings_st';
% path = 'D:\George-abr\human_abr\SingleTrials_aas_04\Emily K. Initial Recordings_st';
% path = 'D:\George-abr\human_abr\SingleTrials_aas_04\Emily K. Repeat Recordings_st';
% path = 'D:\George-abr\human_abr\SingleTrials_aas_05\Emily B. Initial Recordings_st';
path = 'D:\George-abr\human_abr\SingleTrials_aas_05\Emily B. Repeat Recordings_st';

% Obtain directory listing of all ABR files in set (with same electrode
% setup, eg tip horizontal, tip vertical, or wick)
listing = dir(fullfile(path, '*.txt'));
fileList = {listing.name}';
num_files = length(fileList);
is_single_trace = zeros(num_files, 1);
is_montage = zeros(num_files, 1);
is_stack = zeros(num_files, 1);
for i=1:num_files
    filename = fileList{i};
    is_single_trace(i) = contains(filename, SINGLETRACE_KEY, IgnoreCase=true);
    is_montage(i) = contains(filename, MONTAGE, IgnoreCase=true);
    is_stack(i) = is_single_trace(i) && is_montage(i);
end

% Obtain decibel level of files
n_levels = sum(is_stack);
db_levels = zeros(n_levels, 1);
ind_stack = find(is_stack);
is_condensation = zeros(n_levels, 1);
is_rarefaction = zeros(n_levels, 1);
stack = cell(n_levels, 1);
for i=1:n_levels
    filename = fileList{ind_stack(i)};
    ind_underscore = strfind(filename, '_');
    nameNumber = filename(ind_underscore(2) + 1: ind_underscore(3) - 1);
    db_levels(i, 1) = nameToNumber(nameNumber);
    polarity = filename(ind_underscore(3) + 1: ind_underscore(4) - 1);
    is_condensation(i, 1) = strcmpi(polarity, 'C');
    is_rarefaction(i, 1) = strcmpi(polarity, 'R');
    stack{i} = filename;
end

%% Merge single trace files that are of same dB level and opposite polarity
unique_db = unique(db_levels);
for i=1:numel(unique_db)
    this_db = unique_db(i);
    is_db = db_levels==this_db;
    filename_thisdb_condensation = stack{is_db & is_condensation};
    filename_thisdb_rarefaction = stack{is_db & is_rarefaction};
    
    singletrace_cond =  load_humanABR_singletracedata(fullfile(path, filename_thisdb_condensation));
    singletrace_rare =  load_humanABR_singletracedata(fullfile(path, filename_thisdb_rarefaction));
    X_txt = (singletrace_cond + singletrace_rare)/2;
    
    % new filename
    ind_underscore = strfind(filename_thisdb_rarefaction, '_');
    new_filename = [filename_thisdb_rarefaction(1:ind_underscore(3)), 'blk.mat'];
    save(fullfile(path, new_filename), 'X_txt');
end