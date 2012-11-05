function [results,denoiseddata] = GLMdenoisedata(design,data,stimdur,tr,hrfmodel,hrfknobs,opt,figuredir)

% function [results,denoiseddata] = GLMdenoisedata(design,data,stimdur,tr,hrfmodel,hrfknobs,opt,figuredir)
%
% <design> is the experimental design with dimensions time x conditions.
%   Each column should be zeros except for ones indicating condition onsets.
%   Can be a cell vector whose elements correspond to different runs.
%   Different runs can have different numbers of time points.
%   Because this function involves cross-validation across runs, 
%   there must be at least two runs in <design>.
% <data> is the time-series data with dimensions X x Y x Z x time or a cell
%   vector of elements that are each X x Y x Z x time.  XYZ can be collapsed.
%   The dimensions of <data> should mirror that of <design>.  (For example, 
%   <design> and <data> should have the same number of runs, the same number 
%   of time points, etc.)  <data> should not contain any NaNs.
% <stimdur> is the duration of a trial in seconds
% <tr> is the sampling rate in seconds
% <hrfmodel> (optional) indicates the type of model to use for the HRF:
%   'assume' indicates that the HRF is provided (see <hrfknobs>)
%   'optimize' indicates that we should estimate a global HRF from the data
%   Default: 'optimize'.
% <hrfknobs> (optional) is as follows:
%   if <hrfmodel> is 'assume', then <hrfknobs> should be time x 1 with
%     the HRF to assume.
%   if <hrfmodel> is 'optimize', then <hrfknobs> should be time x 1 with the 
%     initial seed for the HRF.  The length of this vector indicates the
%     number of time points that we will attempt to estimate in the HRF.
%   Note on normalization:  We automatically divide <hrfknobs> by the 
%   maximum value so that the peak is equal to 1.  And if <hrfmodel> is 
%   'optimize', then after fitting the HRF, we again normalize the HRF 
%   to peak at 1 (and adjust the amplitudes accordingly).  Default is to 
%   use a canonical HRF that is calculated based on <stimdur> and <tr>.
% <opt> (optional) is a struct with the following fields:
%   <extraregressors> (optional) is time x regressors or a cell vector
%     of elements that are each time x regressors.  The dimensions of 
%     <extraregressors> should mirror that of <design> (i.e. same number of 
%     runs, same number of time points).  The number of extra regressors 
%     does not have to be the same across runs, and each run can have zero 
%     or more extra regressors.  If [] or not supplied, we do 
%     not use extra regressors in the model.
%   <maxpolydeg> (optional) is a non-negative integer with the maximum 
%     polynomial degree to use for polynomial nuisance functions, which
%     are used to capture low-frequency noise fluctuations in each run.
%     Can be a vector with length equal to the number of runs (this
%     allows you to specify different degrees for different runs).  
%     Default is to use round(L/2) for each run where L is the 
%     duration in minutes of a given run.
%   <seed> (optional) is the random number seed to use (this affects
%     the selection of bootstrap samples). Default: sum(100*clock).
%   <bootgroups> (optional) is a vector of positive integers indicating
%     the grouping of runs to use when bootstrapping.  For example, 
%     a grouping of [1 1 1 2 2 2] means that of the six samples that are
%     drawn, three samples will be drawn (with replacement) from the first
%     three runs and three samples will be drawn (with replacement) from
%     the second three runs.  This functionality is useful in situations
%     where different runs involve different conditions.  Default: ones(1,D) 
%     where D is the number of runs.
%   <numforhrf> (optional) is a positive integer indicating the number 
%     of voxels (with the best R^2 values) to consider in fitting the 
%     global HRF.  This input matters only when <hrfmodel> is 'optimize'.
%     Default: 50.  (If there are fewer than that number of voxels
%     available, we just use the voxels that are available.)
%   <brainthresh> (optional) [A B] where A is a percentile for voxel intensity 
%     values and B is a fraction to apply to the percentile.  These parameters
%     are used in the selection of the noise pool.  Default: [99 0.5].
%   <brainR2> (optional) is an R^2 value (percentage).  Voxels whose 
%     cross-validation accuracy is below this value are allowed to enter 
%     the noise pool.  Default: 0.
%   <numpcstotry> (optional) is a non-negative integer indicating the maximum
%     number of PCs to enter into the model.  Default: 10.
%   <numboots> (optional) is a positive integer indicating the number of 
%     bootstraps to perform for the final model.  Special case is 0 which
%     indicates that the final model should just be fit to the complete
%     set of data (and not bootstrapped). Default: 100.
% <figuredir> (optional) is a directory to which to write figures.  (If the
%   directory does not exist, we create it; if the directory already exists,
%   we delete its contents so we can start afresh.)  If [], no figures are
%   written.  If not supplied, default to 'GLMdenoisefigures' (in the current 
%   directory).
%
% Based on the experimental design (<design>, <stimdur>, <tr>) and the 
% model specification (<hrfmodel>, <hrfknobs>), fit a GLM model to the 
% data (<data>, <xyzsize>) using a denoising strategy.  Figures 
% illustrating the results are written to <figuredir>.
%
% Return <results> as a struct containing the following fields:
% <models>, <modelmd>, <modelse>, <R2>, and <R2run> are all like the output of
%   GLMestimatemodel.m (please see that function for details).
% <meanvol> is X x Y x Z with the mean of all volumes
% <noisepool> is X x Y x Z with logical values indicating the voxels that
%   were selected for the noise pool.
% <pcregressors> indicates the global noise regressors that were used
%   to denoise the data.  The format is a cell vector of elements that 
%   are each time x regressors.  The number of regressors will be equal 
%   to opt.numpcstotry.
% <pcR2> is X x Y x Z x (1+opt.numpcstotry) with cross-validated R^2 values for
%   different numbers of PCs.  The first column corresponds to 0 PCs, the
%   second column corresponds to 1 PC, the third column corresponds to
%   2 PCs, etc.
% <pcnum> is the number of PCs that were selected for the final model.
% 
% Also return <denoiseddata>, which is just like <data> except that the 
% component of the data that is estimated to be due to global noise is
% subtracted off.  This may be useful in situations where one wishes to
% treat the denoising as a pre-processing step prior to other analyses 
% of the time-series data.  If you do not need this output, do not
% request it (so we can skip its computation).
%
% Description of the denoising procedure:
% 1. Determine HRF.  If <hrfmodel> is 'assume', we just use the HRF
%    specified by the user.  If <hrfmodel> is 'optimize', we perform
%    a full fit of the GLM model to the data, optimizing the shape of
%    the HRF.
% 2. Determine cross-validated R^2 values.  We fix the HRF to what is
%    obtained in step 1 and estimate the rest of the GLM model.  Leave-one-
%    run-out cross-validation is performed, and we obtain an estimate of the
%    amount of variance (R^2) that can be predicted by the deterministic 
%    portion of the GLM model (the HRF and the amplitudes).
% 3. Determine noise pool.  This is done by calculating a mean volume (the 
%    mean across all volumes) and then determining the voxels that
%    satisfy the following two criteria:
%    (1) The voxels must have sufficient MR signal, that is, the signal
%        intensity in the mean volume must be above a certain threshold
%        (see opt.brainthresh).
%    (2) The voxels must have cross-validated R^2 values that are 
%        below a certain threshold (see opt.brainR2).
% 4. Determine global noise regressors.  This is done by extracting the 
%    time-series data for the voxels in the noise pool, projecting out the
%    polynomial nuisance functions from each time-series, normalizing each
%    time-series to be unit length, and then performing PCA.  The top N
%    PCs from each run (where N is equal to opt.numpcstotry) are selected
%    as global noise regressors.
% 5. Evaluate different numbers of PCs using cross-validation.  We refit
%    the GLM model to the data (keeping the HRF fixed), systematically varying 
%    the number of PCs from 1 to N.  For each number of PCs, leave-one-run-out 
%    cross-validation is performed.  (Recall that only the deterministic
%    portion of the model is cross-validated; thus, any changes in R^2
%    directly reflect changes in the quality of the amplitude estimates.)
% 6. Choose optimal number of PCs.  To choose the optimal number of PCs,
%    we select a simple subset of voxels (namely, all voxels that have 
%    cross-validated R^2 values greater than 0% in the case where no PCs
%    are used) and then compute the median cross-validated R^2 for 
%    these voxels for different numbers of PCs.  The number of PCs with
%    the highest cross-validation accuracy is selected.
% 7. Fit the final model.  We fit the final GLM model (with the HRF fixed to 
%    that obtained in step 1 and with the number of PCs selected in step 6) 
%    to the data.  Bootstrapping is used to estimate error bars on 
%    amplitude estimates.
% 8. Optional: return the de-noised data.  If the output argument <denoiseddata>
%    is requested, we calculate the component of the data that is due to
%    the global noise regressors and return the original time-series data 
%    with this component subtracted off.  Note that the other components of
%    the model (the hemodynamic responses evoked by the experiment, 
%    the polynomial drift, any extra regressors provided by the user) 
%    remain in the de-noised data.
%
% Figures:
% - "HRF.png" shows the initial assumed HRF (provided by <hrfknobs>) and the
%   final estimated HRF (as calculated in step 1).  (If <hrfmodel> is 'assume',
%   the two plots will be identical.)
% - "PCselection.png" shows for different numbers of PCs, the median 
%   cross-validated R^2 across a subset of voxels (namely, those voxels that 
%   have greater than 0% R^2 when no PCs are used).  The optimal number of
%   PCs is indicated in the title of the figure.
% - "PCscatterN.png" shows a scatter plot of cross-validated R^2 values obtained
%   with no PCs against values obtained with N PCs.  The range of the plot
%   is set to the full range of all R^2 values (across all numbers of PCs).
%   To minimize computational time, only up to 50,000 voxels (randomly selected)
%   are shown.
% - "MeanVolume.png" shows the mean volume (mean of all volumes).
% - "NoisePool.png" shows (in white) the voxels selected for the noise pool.
% - "PCcrossvalidationN.png" shows cross-validated R^2 values obtained with N PCs.
%   The range is 0% to 100%, and the colormap is nonlinearly scaled to enhance
%   visibility.
% - "FinalModel.png" shows R^2 values for the final model (as estimated in
%   step 7).  Note that these R^2 values are not cross-validated.
% - "FinalModel_runN.png" shows R^2 values for the final model separated by 
%   runs.  For example, FinalModel_run01.png indicates R^2 values calculated
%   over the data in run 1.
%
% Additional information:
% - For additional details on model estimation and quantification of model
% accuracy (R^2), please see the documentation provided in GLMestimatemodel.m.
%
% History:
% - 2012/11/03 - add a speed-up
% - 2012/11/02 - Initial version.
% - 2012/10/31 - add meanvol and change that it is the mean of all
% - 2012/10/30 - Automatic division of HRF!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEAL WITH INPUTS, ETC.

% input
if ~exist('hrfmodel','var') || isempty(hrfmodel)
  hrfmodel = 'optimize';
end
if ~exist('hrfknobs','var') || isempty(hrfknobs)
  hrfknobs = normalizemax(getcanonicalhrf(stimdur,tr)');
end
if ~exist('opt','var') || isempty(opt)
  opt = struct();
end
if ~exist('figuredir','var')
  figuredir = 'GLMdenoisefigures';
end

% massage input
if ~iscell(design)
  design = {design};
end
if ~iscell(data)
  data = {data};
end

% calc
numruns = length(design);
dataclass = class(data{1});
is3d = size(data{1},4) > 1;
if is3d
  dimdata = 3;
  dimtime = 4;
  xyzsize = sizefull(data{1},3);
else
  dimdata = 1;
  dimtime = 2;
  xyzsize = size(data{1},1);
end

% deal with defaults
if ~isfield(opt,'extraregressors') || isempty(opt.extraregressors)
  opt.extraregressors = cell(1,numruns);
end
if ~isfield(opt,'maxpolydeg') || isempty(opt.maxpolydeg)
  opt.maxpolydeg = zeros(1,numruns);
  for p=1:numruns
    opt.maxpolydeg(p) = round(((size(data{p},dimtime)*tr)/60)/2);
  end
end
if ~isfield(opt,'seed') || isempty(opt.seed)
  opt.seed = sum(100*clock);
end
if ~isfield(opt,'bootgroups') || isempty(opt.bootgroups)
  opt.bootgroups = ones(1,numruns);
end
if ~isfield(opt,'numforhrf') || isempty(opt.numforhrf)
  opt.numforhrf = 50;
end
if ~isfield(opt,'brainthresh') || isempty(opt.brainthresh)
  opt.brainthresh = [99 0.5];
end
if ~isfield(opt,'brainR2') || isempty(opt.brainR2)
  opt.brainR2 = 0;
end
if ~isfield(opt,'numpcstotry') || isempty(opt.numpcstotry)
  opt.numpcstotry = 10;
end
if ~isfield(opt,'numboots') || isempty(opt.numboots)
  opt.numboots = 100;
end
hrfknobs = normalizemax(hrfknobs);
if length(opt.maxpolydeg) == 1
  opt.maxpolydeg = repmat(opt.maxpolydeg,[1 numruns]);
end
if ~iscell(opt.extraregressors)
  opt.extraregressors = {opt.extraregressors};
end

% delete and/or make figuredir
if ~isempty(figuredir)
  if exist(figuredir,'dir')
    assert(rmdir(figuredir,'s'));
  end
  assert(mkdir(figuredir));
  figuredir = absolutepath(figuredir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DETERMINE HRF

% if 'optimize', perform full-fit to determine HRF
if isequal(hrfmodel,'optimize')
  fprintf('*** GLMdenoisedata: performing full fit to estimate global HRF. ***\n');
  fullfit = GLMestimatemodel(design,data,hrfmodel,hrfknobs,0,opt);
  hrf = fullfit.modelmd{1};
  clear fullfit;

% if 'assume', the HRF is provided by the user
else
  hrf = hrfknobs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALCULATE CROSS-VALIDATED R^2 VALUES

% perform cross-validation to determine R^2 values
fprintf('*** GLMdenoisedata: performing cross-validation to determine R^2 values. ***\n');
xvalfit = GLMestimatemodel(design,data,'assume',hrf,-1,opt,1);
pcR2 = xvalfit.R2;
clear xvalfit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DETERMINE NOISE POOL AND CALCULATE GLOBAL NOISE REGRESSORS

% determine noise pool
fprintf('*** GLMdenoisedata: determining noise pool. ***\n');
volcnt = cellfun(@(x) size(x,dimtime),data);
meanvol = reshape(catcell(2,cellfun(@(x) squish(mean(x,dimtime),dimdata),data,'UniformOutput',0)) ...
                  * (volcnt' / sum(volcnt)),[xyzsize 1]);
thresh = prctile(meanvol(:),opt.brainthresh(1))*opt.brainthresh(2);  % threshold for non-brain voxels 
bright = meanvol > thresh;                                           % logical indicating voxels that are bright (brain voxels)
badxval = pcR2 < opt.brainR2;                                        % logical indicating voxels with poor cross-validation accuracy
noisepool = bright & badxval;                                        % logical indicating voxels that satisfy both criteria
  
% determine global noise regressors
fprintf('*** GLMdenoisedata: calculating global noise regressors. ***\n');
pcregressors = {};
for p=1:length(data)

  % extract the time-series data for the noise pool
  temp = subscript(squish(data{p},dimdata),{noisepool ':'})';  % time x voxels

  % project out polynomials from the data
  temp = projectionmatrix(constructpolynomialmatrix(size(temp,1),0:opt.maxpolydeg(p))) * temp;

  % unit-length normalize each time-series (ignoring any time-series that are all 0)
  [temp,len] = unitlengthfast(temp,1);
  temp = temp(:,len~=0);

  % perform SVD and select the top PCs
  [u,s,v] = svds(double(temp*temp'),opt.numpcstotry);
  pcregressors{p} = cast(u,dataclass);

end
clear temp len u s v;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADD GLOBAL NOISE REGRESSORS INTO MODEL AND CHOOSE OPTIMAL NUMBER

% perform cross-validation with increasing number of global noise regressors
for p=1:opt.numpcstotry
  fprintf('*** GLMdenoisedata: performing cross-validation with %d PCs. ***\n',p);
  opt2 = opt;
  for q=1:numruns
    opt2.extraregressors{q} = cat(2,opt2.extraregressors{q},pcregressors{q}(:,1:p));
  end
  xvalfit = GLMestimatemodel(design,data,'assume',hrf,-1,opt2,1);
  pcR2 = cat(dimdata+1,pcR2,xvalfit.R2);
end
clear xvalfit;

% choose optimal number of PCs
temp = squish(pcR2,dimdata);  % voxels x 1+pcs
ok = temp(:,1) > 0;
xvaltrend = median(temp(ok,:),1);
[mx,ix] = max(xvaltrend);
pcnum = ix-1;
fprintf('*** GLMdenoisedata: optimal number of PCs is %d. ***\n',pcnum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIT FINAL MODEL AND PREPARE OUTPUT

% fit final model
opt2 = opt;
for q=1:numruns
  opt2.extraregressors{q} = cat(2,opt2.extraregressors{q},pcregressors{q}(:,1:pcnum));
end
fprintf('*** GLMdenoisedata: fitting final model. ***\n');
results = GLMestimatemodel(design,data,'assume',hrf,opt.numboots,opt2);

% prepare additional outputs
results.meanvol = meanvol;
results.noisepool = noisepool;
results.pcregressors = pcregressors;
results.pcR2 = pcR2;
results.pcnum = pcnum;
clear meanvol noisepool pcregressors pcR2 pcnum;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALCULATE DENOISED DATA

% calculate denoised data (if requested)
if nargout > 1
  fprintf('*** GLMdenoisedata: calculating denoised data. ***\n');
  
  % handle the case of no PCs up front
  if results.pcnum == 0
    denoiseddata = data;

  % otherwise, do regression to figure out the contribution of the PCs
  else 

    % subtract model fit from data
    modelfit = GLMpredictresponses(results.modelmd,design,dimdata);
    data2 = cellfun(@minus,data,modelfit,'UniformOutput',0);
    clear modelfit;
    
    % prepare global noise regressors
    pcmatrix = cellfun(@(x) x(:,1:results.pcnum),results.pcregressors,'UniformOutput',0);
    
    % remove polynomial and other regressors from the global noise regressors and the data
    for p=1:numruns
      temp1 = projectionmatrix(constructpolynomialmatrix(size(data2{p},dimtime),0:opt.maxpolydeg(p)));
      if isempty(opt.extraregressors{p})
        temp2 = 1;
      else
        temp2 = projectionmatrix(opt.extraregressors{p});
      end
      temp = temp1*temp2;
      pcmatrix{p} = temp*pcmatrix{p};
      data2{p} = temp*squish(data2{p},dimdata)';  % time x voxels (flattened format)
    end
    
    % estimate weights on global noise regressors
    h = mtimescell(olsmatrix(blkdiag(pcmatrix{:})),data2);
    clear data2;  % save memory
    
    % construct time-series representing contributions from global noise regressors
    globalnoise = blkdiag(pcmatrix{:})*h;  % time x voxels
    
    % construct denoised data
    denoiseddata = {};
    cnt = 0;
    for p=1:numruns
      ntime = size(data{p},dimtime);
      denoiseddata{p} = data{p} - reshape(globalnoise(cnt+(1:ntime),:)',[xyzsize ntime]);
      cnt = cnt + ntime;
    end
    
    % clean up
    clear pcmatrix globalnoise;

  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERATE FIGURES

if ~isempty(figuredir)
  fprintf('*** GLMdenoisedata: generating figures. ***\n');

  % make figure showing HRF
  figureprep([100 100 450 250]); hold on;
  numinhrf = length(hrfknobs);
  h1 = plot(0:tr:(numinhrf-1)*tr,hrfknobs,'ro-');
  h2 = plot(0:tr:(numinhrf-1)*tr,results.modelmd{1},'bo-');
  ax = axis; axis([0 (numinhrf-1)*tr ax(3) 1.2]);
  straightline(0,'h','k-');
  legend([h1 h2],{'Initial HRF' 'Estimated HRF'});
  xlabel('Time from condition onset (s)');
  ylabel('Response');
  figurewrite('HRF',[],[],figuredir);
  
  % make figure illustrating selection of number of PCs
  figureprep([100 100 400 400]); hold on;
  bar(0:opt.numpcstotry,xvaltrend);
  set(gca,'XTick',0:opt.numpcstotry);
  xlabel('Number of PCs');
  ylabel('Cross-validated R^2 (median across voxels)');
  title(sprintf('Selected PC number = %d',results.pcnum));
  figurewrite('PCselection',[],[],figuredir);
  
  % make figure showing scatter plots of cross-validated R^2
  rng = [min(results.pcR2(:)) max(results.pcR2(:))];
  for p=1:opt.numpcstotry
    temp = squish(results.pcR2,dimdata);  % voxels x 1+pcs
    figureprep([100 100 500 500]); hold on;
    scattersparse(temp(:,1),temp(:,1+p),50000,0,36,'r','.');
    axis([rng rng]); axissquarify; axis([rng rng]); 
    straightline(0,'h','y-');
    straightline(0,'v','y-');
    xlabel('Cross-validated R^2 (0 PCs)');
    ylabel(sprintf('Cross-validated R^2 (%d PCs)',p));
    title(sprintf('Number of PCs = %d (showing at most 50,000 voxels)',p));
    figurewrite(sprintf('PCscatter%02d',p),[],[],figuredir);
  end

  % write out image showing mean volume (of first run)
  imwrite(uint8(255*makeimagestack(results.meanvol,1)),gray(256),[figuredir '/MeanVolume.png']);

  % write out image showing noise pool
  imwrite(uint8(255*makeimagestack(results.noisepool,[0 1])),gray(256),[figuredir '/NoisePool.png']);

  % define a function that will write out R^2 values to an image file
  imfun = @(results,filename) ...
    imwrite(uint8(255*makeimagestack(signedarraypower(results/100,0.5),[0 1])),hot(256),filename);

  % write out cross-validated R^2 for the various numbers of PCs
  for p=1:size(results.pcR2,dimdata+1)
    temp = subscript(results.pcR2,[repmat({':'},[1 dimdata]) {p}]);
    feval(imfun,temp,sprintf([figuredir '/PCcrossvalidation%02d.png'],p-1));
  end

  % write out overall R^2 for final model
  feval(imfun,results.R2,sprintf([figuredir '/FinalModel.png']));

  % write out R^2 separated by runs for final model
  for p=1:size(results.R2run,dimdata+1)
    temp = subscript(results.R2run,[repmat({':'},[1 dimdata]) {p}]);
    feval(imfun,temp,sprintf([figuredir '/FinalModel_run%02d.png'],p));
  end
  
end
