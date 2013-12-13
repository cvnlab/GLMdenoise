%% Example 3: Run GLMdenoise using finite impulse response (FIR) model

%% Download dataset (if necessary) and add GLMdenoise to the MATLAB path

setup;

%% Load in the data

% Load in the data
load('exampledataset.mat');

% Check the workspace
whos
%%

%% Call GLMdenoisedata using FIR model

% It is typical to assume that the same fixed hemodynamic response function (HRF) 
% applies across all experimental conditions and all voxels.  Indeed, this is the 
% default behavior of GLMdenoisedata.
%
% In contrast, the finite impulse response (FIR) model makes no assumption about the
% shape of the hemodynamic response; it does this by estimating a separate parameter 
% for each time point in the response to each condition.
%
% Applying the FIR model can be useful for datasets where differences in timecourses
% exist across conditions or voxels.  Moreover, even if timecourse differences do not
% necessarily exist, it may still be useful to apply the FIR model as a way of checking
% that this is indeed the case.  Also, applying the FIR model can be used as a means
% of troubleshooting (e.g., it may reveal some error in the handling or preparation 
% of the data).
%
% GLMdenoisedata includes a flag that allows the FIR model to be used.  Denoising
% proceeds as usual because the FIR model simply changes the model of the experiment-
% related effects in the data.

% Here we apply the FIR model to the example dataset.  The 'fir' option indicates to
% use the FIR model and the 20 indicates to estimate timecourses up to 20 time 
% points after condition onset.  Notice that the design matrix passed to GLMdenoisedata 
% is in the usual format (a matrix that indicates condition onsets); there is no 
% need to construct the FIR basis functions as that is done internally by the code.
% Also, notice that here we turn off bootstrapping to save computational time and 
% memory requirements.
results = GLMdenoisedata(design,data,stimdur,tr, ...
            'fir',20,struct('numboots',0), ...
            'example3figures');
%%

% Let's inspect the results for an example voxel.  By eye, it appears that the 
% FIR timecourses for this voxel are reasonably approximated by a single HRF.
% Also, notice that the canonical HRF (red line) is substantially different from 
% the FIR timecourses, suggesting that there may be value in fitting an HRF to the data.
xx = 53; yy = 27; zz = 1;
tcs = squeeze(results.modelmd(xx,yy,zz,:,:));  % 35 conditions x 21 time points
figure; hold on;
plot(0:tr:20*tr,tcs');  % plot FIR timecourses
plot(0:tr:20*tr,mean(tcs,1),'k-','LineWidth',3);  % plot the mean timecourse
hrf = getcanonicalhrf(stimdur,tr);  % compute a canonical HRF
  % scale canonical HRF to best match the mean timecourse
hrf = hrf * (pinv(hrf(1:21)')*mean(tcs,1)');
plot(0:tr:(length(hrf)-1)*tr,hrf,'r-','LineWidth',3);
straightline(0,'h','k-');
ax = axis;
axis([0 30 ax(3:4)]);
xlabel('Time from condition onset (s)');
ylabel('BOLD signal (% change)');
title('FIR timecourses (black = mean, red = canonical HRF)');
%%
