function [data, session_name] = bids_loadfscans(sub,prep_path,session_type, session_number, structural_space)
% 2020-06-04 written by Alex Kipnis (tested on https://openneuro.org/datasets/ds001246/versions/1.2.1 after running fMRIprep on sub-01) 
%     Load functional scans (nii.gz files) into "data" cell-array, used as
%     input for GLMdenoise.
%
%   
%     Parameters: 
%     sub (int): subject number
%     prep_path (str): path to preprocessed data (typically in
%      "derivatives" subfolder 
%     session_type (str): all folders containing this substring will be
%       searched for nii.gz files
%     session_number (int): Self-explanatory (GLMs should be fit with
%       sessions and runs as separate levels)
%     structural_space (str): e.g. T1w or MNI152
%   
%     Returns: 
%     data (cell-array): Nifti 4D-arrays for each run
%     session_name (str): title for the session from which the runs are imported

% Initialize variables
n = session_number; file_list = []; filename_list = {}; data = {};
fprintf('Searching for subject %d''s NIfTI files, containing functional scans for "%s" sessions in %s space...\n',...
    sub, session_type, structural_space)

% Create path structures
type_path = strcat(prep_path, string(sub), filesep ,'ses-', session_type,'*');
dir_prep = dir(type_path); %folder list of preprocessed data (we need this for nifti files);
dir_ses = dir(strcat(dir_prep(n).folder, filesep, dir_prep(n).name, filesep, 'func',filesep)); %folder list for specified number of the session type
session_name = dir_prep(n).name;

% collect file names for individual runs 
for f = 1:length(dir_ses)
   if contains(dir_ses(f).name, structural_space) & contains(dir_ses(f).name,'bold.nii.gz')
       file_list = [file_list; f]; %append list when filename matches conditions
       filename_list = cat(2, filename_list, dir_ses(f).name);
   end
end

% load files from filename_list
for g = 1:length(file_list)
    fprintf('Loading run %d... \n', g)
    run = niftiread(strcat(dir_ses(f).folder, '/', string(filename_list(g)))); %append folder name with file name and load nii.gz file
    data = cat(2, data, run);
end

