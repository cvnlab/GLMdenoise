function [design, stimdur, ints2unique, unique2ints] = bids_loaddesign(sub,BIDS_path,session_type, session_number, tr)
%%% 2020-05-30 written by Alex Kipnis (tested on https://openneuro.org/datasets/ds001246/versions/1.2.1 after running fMRIprep on sub-01) 
%     Loads event markers (tsv files) into "design" cell-array, which is used as
%     input for GLMdenoise. This script assumes that the presentation
%     duration for each stimulus is constant and regressors that are not
%     stimulus IDs (e.g. prespecified noise regressors) will be added to the design matrix separately. 
%
%   
%     Parameters: 
%     sub (int): subject number
%     BIDS_path (str): path to BIDS-formatted folder 
%     session_type (str): all folders containing this substring will be
%       searched for .tsv files
%     session_number (int): Self-explanatory (GLMs should be fit with
%       sessions and runs as separate levels)
%     tr (int): Repetition Time in seconds
%   
%     Returns: 
%     design (sparse-double): sparse samples x stim-type matrices
%     stimdur (int): stimulus duration in seconds
%     ints2unique(dictionary): map between keys (integers of individual
%       stimuli) and values (stimulus names)
%     unique2ints(dictionary): self-explanatory

% Initialize variables
n = session_number; file_list = []; filename_list = {}; all_events={}; design = {}; stim_ids = char(); ints2unique={};
fprintf('Transforming subject %d''s .tsv event files for "%s" sessions into design matrices.\n', sub, session_type)

% Create path structures
type_path = strcat(BIDS_path, 'sub-0', string(sub), '/ses-', session_type,'*');
dir_BIDS = dir(type_path); %folder list of preprocessed data (we need this for nifti files);
dir_ses = dir(strcat(dir_BIDS(n).folder, '/', dir_BIDS(n).name, '/func/')); %folder list for specified number of the session type

% collect file names for individual runs 
for f = 1:length(dir_ses)
   if contains(dir_ses(f).name,'events.tsv')
       file_list = [file_list; f]; %append list when filename matches conditions
       filename_list = cat(2, filename_list, dir_ses(f).name);
   end
end

% load files from filename_list
for g = 1:length(file_list)
    events = tdfread(strcat(dir_ses(f).folder, '/', string(filename_list(g)))); %append folder name with file name and load .tsv file    
    all_events = cat(2, all_events, events); %save all events
    stim_ids = cat(1, stim_ids, events.stim_id); %concatenate all stimulus IDs
end

% save stimulus duration
stimdur = events.duration(2); %take the second duration (as the first is always just the beginning of the run);

% make list of unique events
unique_ids = unique(cellstr(stim_ids));
idx=cellfun(@(x) isequal(x,'n/a'),unique_ids); %indices of 'n/a' entries
unique_ids(idx,:)=[]; %remove it.
integer_ids = 1:length(unique_ids);
ints2unique = containers.Map(integer_ids, unique_ids); % create event dictionary
unique2ints = containers.Map(unique_ids, integer_ids); % create event dictionary

% make sparse matrices
if stimdur/tr <= 1 %if the stimulus duration was not longer than time between two consecutive volume samples
    for e = 1:length(all_events)
    %     trial_no = all_events{e}.trial_no; trial_no(1) = []; trial_no(end) = []; % extract trial numbers
        stim_onset = all_events{e}.onset; stim_onset(1,:) = []; stim_onset(end,:) = []; % extract stimulus onset
        stim_sample = stim_onset/tr; % transform onset timepoint to sampling point by dividing by TR
        stim_id = cellstr(all_events{e}.stim_id); stim_id(1,:) = []; stim_id(end,:) = []; % extract stimulus IDs
        event_vector = values(unique2ints,stim_id); % map stimulus IDs onto their corresponding integers
        event_vector = cell2mat(event_vector);
        M = sparse(stim_sample, event_vector, ones(length(event_vector),1)); %create sparse matrix with with ones for each tuple of event and sampling point
        design{e} = M; %save all design matrices
    end
else %if a stimulus was presented over several volume samples 
    for e = 1:length(all_events)
        sampleperstim = stimdur/tr;
    %     trial_no = all_events{e}.trial_no; trial_no(1) = []; trial_no(end) = []; % extract trial numbers
        stim_onset = all_events{e}.onset; stim_onset(1,:) = []; stim_onset(end,:) = []; % extract stimulus onset
        stim_sample = stim_onset/tr; % transform onset timepoint to sampling point by dividing by TR
        stim_id = cellstr(all_events{e}.stim_id); stim_id(1,:) = []; stim_id(end,:) = []; % extract stimulus IDs
        event_vector = values(unique2ints,stim_id); % map stimulus IDs onto their corresponding integers
        event_vector = cell2mat(event_vector);
       
        
        %account for the fact that each stimulus is presented over several sampled volumes
        stim_sample_orig = stim_sample; event_vector_orig = event_vector;
        for s = 1:(sampleperstim-1)
            stim_sample_temp = stim_sample_orig+s; %add s to each samplepoint 
            stim_sample = cat(1, stim_sample, stim_sample_temp); %add the new samples to previous sample_point vector
            event_vector = cat(1, event_vector, event_vector_orig); %"multiply" event_vector accordingly
        end   
        %create sparse matrix (of size m total sampling points x n distinct stimuli) with with ones for each tuple of stimulus and sampling point
        M = sparse(stim_sample, event_vector, ones(length(event_vector),1), sum(events.duration)/tr, length(unique_ids)); 
        design{e} = M; %save all design matrices
    end
end
