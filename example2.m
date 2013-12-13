%% Example 2: Analyze different splits of runs in a session

%% Download dataset (if necessary) and add GLMdenoise to the MATLAB path

setup;

%% Load in the data

% Load in the data
load('exampledataset.mat');

% Check the workspace
whos
%%

%% Split session into two halves and analyze each half separately

% In some circumstances, you may want to split the runs acquired in a given
% session into groups and analyze each group separately.  For example, you might
% want to train a classifier on data from some runs and then test the 
% classifer on data from other runs.  Here we go through an example where
% we split a session of data (10 runs) into halves and analyze each half separately.

% First, we define the way the runs are to be split.  In the first
% iteration, we will analyze the odd runs.  In the second iteration,
% we will analyze the even runs.
splits = {[1:2:10] [2:2:10]};

% For each split, we will analyze the data two ways.  First, we analyze the data
% using no denoising (REGULAR).  Then, we analyze the data using denoising (DENOISE).
% (To achieve no denoising, we set the number of PCs to try to zero.)
clear resultsREGULAR resultsDENOISE;
for p=1:length(splits)
  ix = splits{p};

  % analyze using no denoising
  resultsREGULAR(p) = GLMdenoisedata(design(ix),data(ix),stimdur,tr, ...
                        [],[],struct('numpcstotry',0), ...
                        sprintf('example2figures_REGULAR_split%d',p));

  % analyze using denoising
  resultsDENOISE(p) = GLMdenoisedata(design(ix),data(ix),stimdur,tr, ...
                        [],[],[], ...
                        sprintf('example2figures_DENOISE_split%d',p));

end
%%

% The default mode of operation involves estimating an HRF from the data (in each
% call to GLMdenoisedata).  Let's inspect the HRF results across the two splits of 
% the data.  The reason this may be important is that the beta weights that are 
% estimated in the model are scale factors that are applied to the HRF estimate,
% and if the HRF estimates are highly different across the splits, that may be some
% reason for concern.  In this case, we find that the HRF estimate is fairly 
% similar across the splits.  (An alternative strategy is to just fix the HRF to
% a canonical HRF and use that for each split.)
figure; hold on;
hrf1 = resultsREGULAR(1).modelmd{1};
hrf2 = resultsREGULAR(2).modelmd{1};
plot(0:tr:tr*(length(hrf1)-1),hrf1,'r-');
plot(0:tr:tr*(length(hrf2)-1),hrf2,'b-');
xlabel('Time from condition onset (s)');
ylabel('Response');
title('HRF estimate across splits (red = odd runs, blue = even runs)');
%%

% Let's now examine the stability of beta weights across the data splits.
% First, we choose a subset of the voxels.
ix = find(resultsREGULAR(1).R2 > 10 & resultsREGULAR(1).meanvol > 500);

% Then, we extract the beta weights corresponding to one of the
% conditions (condition 10).  Notice that we have four sets of beta weights.
% This is because for each of the two splits of the data, we have two 
% beta weight estimates (one using no denoising, one using denoising). 
beta1REGULAR = subscript(squish(resultsREGULAR(1).modelmd{2},3),{ix 10});
beta2REGULAR = subscript(squish(resultsREGULAR(2).modelmd{2},3),{ix 10});
beta1DENOISE = subscript(squish(resultsDENOISE(1).modelmd{2},3),{ix 10});
beta2DENOISE = subscript(squish(resultsDENOISE(2).modelmd{2},3),{ix 10});

% Visualize the results.
figure;
set(gcf,'Units','points','Position',[100 100 800 400]);
subplot(1,2,1); hold on;
scatter(beta1REGULAR,beta2REGULAR,'r.');
axissquarify;
xlabel('BOLD signal (% change), odd runs');
ylabel('BOLD signal (% change), even runs');
title(sprintf('No denoising (r=%.3f)',corr(beta1REGULAR,beta2REGULAR)));
subplot(1,2,2); hold on;
scatter(beta1DENOISE,beta2DENOISE,'r.');
axissquarify;
xlabel('BOLD signal (% change), odd runs');
ylabel('BOLD signal (% change), even runs');
title(sprintf('Denoising (r=%.3f)',corr(beta1DENOISE,beta2DENOISE)));
%%

%% Split session into individual runs and analyze each run separately

% In some circumstances, you may wish to derive beta weights from individual
% runs in a session.  Here we go through an example where we analyze each
% of the runs in a session separately.  Note that denoising cannot be applied
% if there is only a single run to work with (since cross-validation is performed
% across runs).  So, instead of calling GLMdenoisedata.m, we will call
% GLMestimatemodel.m.  Calling GLMestimatemodel.m also has the benefit that 
% we will avoid some of the overhead involved in GLMdenoisedata.m.

% First, we define the ways the runs are to be split.  Every run will be
% analyzed separately.
splits = num2cell(1:10);

% For demonstration purposes, we will adopt the strategy of assuming a fixed HRF
% in the analysis of each run.  Here we compute a canonical HRF.
hrf = getcanonicalhrf(stimdur,tr)';

% Now we loop over each run, using GLMestimatemodel.m to fit a GLM to the data.
% Notice that we specify that we want to use an assumed HRF.  Also, the "0"
% indicates that we just want to fit the data (no bootstrapping nor cross-validation).
clear resultsIND;
for p=1:length(splits)
  ix = splits{p};
  resultsIND(p) = GLMestimatemodel(design(ix),data(ix),stimdur,tr,'assume',hrf,0);
end
%%

% Let's examine the beta weights from each run for an example voxel. 
xx = 46; yy = 24; zz = 1;
figure; hold on;
set(gcf,'Units','points','Position',[100 100 700 250]);
betas = [];
for p=1:length(resultsIND)
  betas(p,:) = flatten(resultsIND(p).modelmd{2}(xx,yy,zz,:));
  plot(betas(p,:),'b-');
end
plot(mean(betas,1),'k-','LineWidth',3);
straightline(0,'h','k-');
xlabel('Condition number');
ylabel('BOLD signal (% change)');
title('Beta weights (blue = estimates from individual runs, black = mean across estimates)');
%%
