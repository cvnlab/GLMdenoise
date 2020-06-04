%% Example 5: Run GLMdenoise on data in BIDS format
%% 2020-05-30 written by Alex Kipnis (tested on https://openneuro.org/datasets/ds001246/versions/1.2.1 after running fMRIprep on sub-01)
% This script assumes you have a raw dataset in BIDS format (https://bids.neuroimaging.io/) with preprocessed data in its "derivatives" subdirectory 

% 0. Specify Path variables
GLMd_path = "/home/alex/Documents/MATLAB/GLMdenoise"; cd(GLMd_path); %path to GLMdenoise toolbox
BIDS_path = "/home/alex/ds001246/"; %folder path to BIDS formatted data
prep_path = strcat('/mnt/locker/alex/ds001246/derivatives/fmriprep/sub-0'); %subfolder path to preprocessed data
json_path = strcat(BIDS_path, 'task-perception_bold.json'); %path to metadata (.json file)

%% 1. Load in the data
sub = 1; %subject id
session_type = 'perception';
session_number = 2;
structural_space = "T1w";
metadata = jsondecode(fileread(json_path)); tr = metadata.RepetitionTime;

% 1.1 Load nifti files
[data, session_name] = bids_loadfscans(sub,prep_path,session_type, session_number, structural_space);

% 1.2 Load tsv files
[design, stimdur, event_dict] = bids_loaddesign(sub,BIDS_path,session_type, session_number, tr);

%% Inspect the data

% Look at the design matrix for the first run (50 possible images, one
% being presented at the first sampling point of each trial)
figure;
set(gcf,'Units','points','Position',[100 100 500 1770]);
imagesc(design{1});
colormap(gray);
colorbar;
xlabel('Image');
ylabel('Time points');
title(sprintf('Design matrix (%d sampling points, %d conditions)',size(design{1},1),size(design{1},2)));

%%
% Look at the first volume of the first run.  We split the volume along
% the third dimension (the slice dimension) and show slices starting at the
% upper-left, proceeding downwards, and then from left to right.
figure;
imagesc(makeimagestack(data{1}(:,:,:,1)));
colormap(gray);
axis equal tight off;
colorbar;
title('fMRI data (first volume)');
%%

% Check various constants
fprintf('There are %d runs in total.\n',length(design));
fprintf('The dimensions of the data for the first run are %s.\n',mat2str(size(data{1})));
fprintf('The stimulus duration is %.6f seconds.\n',stimdur);
fprintf('The sampling rate (TR) is %.6f seconds.\n',tr);
%%


%% Run GLMdenoise

% Call GLMdenoise using default parameters.  Outputs will be stored in
% 'results' and 'denoiseddata'.  Figures will be written to 'example5figures'
% in the current directory.
[results,denoiseddata] = GLMdenoisedata(design,data,stimdur,tr,[],[],struct('numpcstotry',1),'example5figures');



% Sanity check
 
%% Inspect figures
% 'FinalModel.png' shows the R^2 of the final model fit.  The color range is
% 0% to 100% (but is nonlinear; see the color bar).  These R^2 values are not 
% cross-validated (thus, even voxels with no signal have positive R^2 values).
figure;
imagesc(imread('kamitanifigures/FinalModel.png'),[0 255]);
colormap(hot);
axis equal tight off;
cb = colorbar;
set(cb,'YTick',linspace(0,255,11),'YTickLabel',round((0:.1:1).^2 * 100));
title('Final model R^2 (not cross-validated)');
%%

% The outputs of GLMdenoisedata are contained in 
% the variables 'results' and 'denoiseddata'.
% Here we do some basic inspections of the outputs.

% Select a voxel to inspect.  This is done by finding voxels that have 
% cross-validated R^2 values between 0% and 5% under the initial model (no PCs),
% and then selecting the voxel that shows the largest improvement when
% using the final model.
ix = find(results.pcR2(:,:,:,1) > 0 & results.pcR2(:,:,:,1) < 5);
improvement = results.pcR2(:,:,:,1+results.pcnum) - results.pcR2(:,:,:,1);
[mm,ii] = max(improvement(ix));
ix = ix(ii);
[xx,yy,zz] = ind2sub(results.inputs.datasize{1}(1:3),ix);


% An alternative to using the GLM estimates provided by GLMdenoisedata is
% to use 'denoiseddata', which contains the original time-series data but 
% with the component of the data that is estimated to be due to the global 
% noise regressors subtracted off.  Here we inspect the denoised data for 
% the same example voxel examined earlier.  Note that the data components 
% that are present in 'denoiseddata' can be customized (see opt.denoisespec 
% in GLMdenoisedata.m).  The default parameters leave in the estimated baseline
% signal drift, which explains the drift in the plotted time-series.
figure; hold on;
set(gcf,'Units','points','Position',[100 100 700 250]);
data1 = flatten(data{1}(xx,yy,zz,:));
data2 = flatten(denoiseddata{1}(xx,yy,zz,:));
n = length(data1);
h1 = plot(data1,'r-');
h2 = plot(data2,'b-');
ax = axis; axis([0 n+1 ax(3:4)]);
legend([h1 h2],{'Original' 'Denoised'});
xlabel('Time point');
ylabel('MR signal');

% Optionally: Save data and clear your workspace
s = struct('subject_id',[sub], 'session_type', [session_type], 'session_number', [session_number],...
    'session_name', [session_name], 'structural_space', [structural_space], 'metadata', [metadata],...
    'n_runs', [length(design)], 'data_dims', [mat2str(size(data{1}))], 'stimdur', [stimdur],...
    'design', [{design}], 'event_dict', [event_dict], 'mean_volume', [results.meanvol],...
    'HRF', [results.modelmd{1, 1}], 'beta_coefficients', [results.modelmd{1, 2}],...
    'cross_validated_R2', [results.pcR2], 'number_of_principle_components',[results.pcnum]);

clearvars -except s
names = fieldnames(s);
for i = 1:numel(names)
    assignin('caller', names{i}, s.(names{i}));
end
clear s names i

filename = strcat('sub-0', string(subject_id), '-', session_name,'-export.mat');
save(filename)
