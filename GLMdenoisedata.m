function [results,denoiseddata] = GLMdenoisedata(design,data,stimdur,tr,hrfmodel,hrfknobs,opt,figuredir)

% function [results,denoiseddata] = GLMdenoisedata(design,data,stimdur,tr,hrfmodel,hrfknobs,opt,figuredir)
%
% <design> is the experimental design.  There are three possible cases:
%   1. A where A is a matrix with dimensions time x conditions.
%      Each column should be zeros except for ones indicating condition onsets.
%      (Fractional values in the design matrix are also allowed.)
%   2. {A1 A2 A3 ...} where each of the A's are like the previous case.
%      The different A's correspond to different runs, and different runs
%      can have different numbers of time points.
%   3. {{C1_1 C2_1 C3_1 ...} {C1_2 C2_2 C3_2 ...} ...} where Ca_b is a vector of 
%      onset times (in seconds) for condition a in run b.  Time starts at 0 
%      and is coincident with the acquisition of the first volume.  This case 
%      is compatible only with <hrfmodel> set to 'assume'.
%   Because this function involves cross-validation across runs, there must 
%   be at least two runs in <design>.  
% <data> is the time-series data with dimensions X x Y x Z x time or a cell vector of 
%   elements that are each X x Y x Z x time.  XYZ can be collapsed such that the data 
%   are given as a 2D matrix (XYZ x time); however, if you do this, then several of the 
%   figures that are written out by this function will not be useful to look at.  The 
%   dimensions of <data> should mirror that of <design>.  (For example, <design> and 
%   <data> should have the same number of runs, the same number of time points, etc.  
%   Thus, <data> should have at least two runs since <design> must have at least two 
%   runs.)  <data> should not contain any NaNs. We automatically convert <data> to 
%   single format (if necessary).
% <stimdur> is the duration of a trial in seconds
% <tr> is the sampling rate in seconds
% <hrfmodel> (optional) indicates the type of model to use for the HRF:
%   'fir' indicates a finite impulse response model (a separate timecourse
%     is estimated for every voxel and every condition)
%   'assume' indicates that the HRF is provided (see <hrfknobs>)
%   'optimize' indicates that we should estimate a global HRF from the data
%   Default: 'optimize'.
% <hrfknobs> (optional) is as follows:
%   if <hrfmodel> is 'fir', then <hrfknobs> should be the number of 
%     time points in the future to model (N >= 0).  For example, if N is 10, 
%     then timecourses will consist of 11 points, with the first point 
%     coinciding with condition onset.
%   if <hrfmodel> is 'assume', then <hrfknobs> should be time x 1 with
%     the HRF to assume.
%   if <hrfmodel> is 'optimize', then <hrfknobs> should be time x 1 with the 
%     initial seed for the HRF.  The length of this vector indicates the
%     number of time points that we will attempt to estimate in the HRF.
%   Note on normalization:  In the case that <hrfmodel> is 'assume' or
%   'optimize', we automatically divide <hrfknobs> by the maximum value
%   so that the peak is equal to 1.  And if <hrfmodel> is 'optimize',
%   then after fitting the HRF, we again normalize the HRF to peak at 1
%   (and adjust amplitudes accordingly).  Default in the case of 'fir' is
%   20.  Default in the case of 'assume' and 'optimize' is to use a 
%   canonical HRF that is calculated based on <stimdur> and <tr>.
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
%   <hrffitmask> (optional) is X x Y x Z with 1s indicating all possible
%     voxels to consider for fitting the global HRF.  This input matters
%     only when <hrfmodel> is 'optimize'.  Special case is 1 which means
%     all voxels can be potentially chosen.  Default: 1.
%   <brainthresh> (optional) [A B] where A is a percentile for voxel intensity 
%     values and B is a fraction to apply to the percentile.  These parameters
%     are used in the selection of the noise pool.  Default: [99 0.5].
%   <brainR2> (optional) is an R^2 value (percentage).  Voxels whose 
%     cross-validation accuracy is below this value are allowed to enter 
%     the noise pool.  Default: 0.
%   <brainexclude> (optional) is X x Y x Z with 1s indicating voxels to 
%     exclude when selecting the noise pool.  Special case is 0 which means 
%     all voxels can be potentially chosen.  Default: 0.
%   <numpcstotry> (optional) is a non-negative integer indicating the maximum
%     number of PCs to enter into the model.  Default: 20.
%   <pcR2cutoff> (optional) is an R^2 value (percentage).  To decide the number
%     of PCs to include, we examine a subset of the available voxels.
%     Specifically, we examine voxels whose cross-validation accuracy is above 
%     <pcR2cutoff> for any of the numbers of PCs.  Default: 0.
%   <pcR2cutoffmask> (optional) is X x Y x Z with 1s indicating all possible
%     voxels to consider when selecting the subset of voxels.  Special case is
%     1 which means all voxels can be potentially selected.  Default: 1.
%   <pcstop> (optional) is
%     A: a number greater than or equal to 1 indicating when to stop adding PCs 
%        into the model.  For example, 1.05 means that if the cross-validation 
%        performance with the current number of PCs is within 5% of the maximum 
%        observed, then use that number of PCs.  (Performance is measured 
%        relative to the case of 0 PCs.)  When <pcstop> is 1, the selection 
%        strategy reduces to simply choosing the PC number that achieves
%        the maximum.  The advantage of stopping early is to achieve a selection
%        strategy that is robust to noise and shallow performance curves and 
%        that avoids overfitting.
%    -B: where B is the number of PCs to use for the final model (thus, the user
%        chooses).  B can be any integer between 0 and opt.numpcstotry.
%     Default: 1.05.
%   <pccontrolmode> (optional) is for testing purposes.  Default is 0 which means
%     to do nothing special.  If 1, then after we are done constructing the global
%     noise regressors, we scramble the phase spectra of these regressors (prior
%     to entering them into the GLM).  If 2, then after we are done constructing the
%     noise regressors, we shuffle the assignment of noise regressors
%     to the runs, ensuring that each run is assigned a new set of regressors.  Note that
%     in this shuffling, the grouping of regressors (at the run level) is maintained.
%     The shuffling is performed prior to entering noise regressors into the GLM.
%   <numboots> (optional) is a positive integer indicating the number of 
%     bootstraps to perform for the final model.  Special case is 0 which
%     indicates that the final model should just be fit to the complete
%     set of data (and not bootstrapped). Default: 100.
%   <denoisespec> (optional) is a binary string or cell vector of binary strings
%     indicating the components of the data to return in <denoiseddata>.  The 
%     format of each string should be 'ABCDE' where A indicates whether to include 
%     the signal (estimated hemodynamic responses evoked by the experiment), B 
%     indicates whether to include the polynomial drift, C indicates whether
%     to include any extra regressors provided by the user, D indicates 
%     whether to include the noise regressors, and E indicates whether
%     to include the residuals of the model.  If multiple strings are provided,
%     then separate copies of the data will be returned in the rows of 
%     <denoiseddata>.  Default: '11101' which indicates that all components of 
%     the data will be returned except for the component corresponding to the 
%     estimate of the contribution of the noise regressors.
%   <wantpercentbold> (optional) is whether to convert the amplitude estimates
%     in 'models', 'modelmd', and 'modelse' to percent BOLD change.  This is
%     done as the very last step, and is accomplished by dividing by the 
%     absolute value of 'meanvol' and multiplying by 100.  (The absolute 
%     value prevents negative values in 'meanvol' from flipping the sign.)
%     Default: 1.
%   <hrfthresh> (optional) is an R^2 threshold.  If the R^2 between the estimated
%     HRF and the initial HRF is less than <hrfthresh>, we decide to just use
%     the initial HRF.  Set <hrfthresh> to -Inf if you never want to reject
%     the estimated HRF.  Default: 50.
%   <noisepooldirect> (optional) is {A B} where A is X x Y x Z with 1s indicating the 
%     voxels to use as the noise pool and B is a non-negative integer indicating the 
%     number of PCs to use.  (B must be less than or equal to <numpcstotry>.)
%     If <noisepooldirect> is supplied, this causes much of the standard GLMdenoise 
%     procedure to be bypassed.  For example, cross-validation is no longer necessary 
%     and therefore no longer performed.  Thus, one benefit of using <noisepooldirect>
%     is that you can apply GLMdenoisedata.m to a single run of data.  Note that if 
%     <noisepooldirect> is supplied, various inputs no longer have any effect 
%     (e.g., <brainthresh>, <brainR2>, <brainexclude>, <pcR2cutoff>, <pcR2cutoffmask>, 
%     and <pcstop>) and various outputs and figures are omittted.  Default is [] which 
%     means to perform the usual GLMdenoise procedure.
%   <wantparametric> (optional) is whether to compute parametric GLM fits and associated
%     error estimates.  These are added as additional outputs in the <results> struct.
%     Default: 0.
%   <wantsanityfigures> (optional) is whether to write out figures that allow you
%     to sanity-check the data.  You may want to set this to 0 to save computational
%     time.  Default: 1.
%	  <drawfunction> (optional) is a function that accepts values that are X x Y x Z
%     (or XYZ x 1) and returns a 2D or 3D matrix. This matrix is then passed to
%     makeimagestack.m for the purposes of writing image files. Specifying 
%     <drawfunction> can be useful for transforming values that are
%     column vectors (XYZ x 1) into a more palatable form. Default is to
%     do nothing special: @(vals) vals. Note that we do not save this input
%     to results.inputs.opt in order to avoid weird .mat file-saving problems.
%   <reconmask> (optional) is a X x Y x Z mask that can be used to reconstruct 3D 
%     images just prior to figure creation. This can only be done if the data
%     is provided in vector (XYZ x Time) format AND if the vector doesn't contain 
%     any zeros (a mask was used to create it). This 'reduced vector' only contains 
%     data, and has no zero timeseries in it - massively reduces RAM required.
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
% <models>, <modelmd>, <modelse>, <R2>, <R2run>, <signal>, <noise>, 
%   <SNR>, and <hrffitvoxels> are all like the output of GLMestimatemodel.m 
%   (please see that function for details).
% <hrffitvoxels> is X x Y x Z with logical values indicating the voxels that
%   were used to fit the global HRF.  (If <hrfmodel> is not 'optimize',
%   this is returned as [].)
% <meanvol> is X x Y x Z with the mean of all volumes
% <noisepool> is X x Y x Z with logical values indicating the voxels that
%   were selected for the noise pool.
% <pcregressors> indicates the noise regressors that were used
%   to denoise the data.  The format is a cell vector of elements that 
%   are each time x regressors.  The number of regressors will be equal 
%   to opt.numpcstotry and they are in order (the first PC is the first
%   regressor, the second PC is the second regressor, etc.).  Note that
%   only the first <pcnum> PCs are actually used for the final model.
% <pcR2> is X x Y x Z x (1+opt.numpcstotry) with cross-validated R^2 values for
%   different numbers of PCs.  The first slice corresponds to 0 PCs, the
%   second slice corresponds to 1 PC, the third slice corresponds to
%   2 PCs, etc.
% <pcR2final> is X x Y x Z with cross-validated R^2 values for the final
%   model that is chosen.  Note that this is simply pcR2(:,:,:,1+pcnum).
% <pcvoxels> is X x Y x Z with logical values indicating the voxels that
%   were used to select the number of PCs.
% <pcnum> is the number of PCs that were selected for the final model.
% <pcweights> is X x Y x Z x <pcnum> x R with the estimated weights on the 
%   PCs for each voxel and run.
% <SNRbefore> is X x Y x Z with signal-to-noise ratios before denoising
%   (using no noise regressors).
% <SNRafter> is X x Y x Z with signal-to-noise ratios after denoising
%   (using the final number of noise regressors).
% <inputs> is a struct containing all inputs used in the call to this
%   function, excluding <data>.  We additionally include a field called 
%   'datasize' which contains the size of each element of <data>.
% <parametric> is a struct that is included if opt.wantparametric.  The fields are:
%   <designmatrix> is time x regressors with the full design matrix of the GLM
%   <parameters> is X x Y x Z x regressors with the beta weight estimate for
%     each regressor
%   <parametersse> is X x Y x Z x regressors with the standard error on the beta weights
%   <noisevar> is X x Y x Z with the estimate of the noise variance (sigma squared)
%   Note that these are in raw units and are not converted into % BOLD change.
% 
% Also return <denoiseddata>, which is just like <data> except that the 
% component of the data that is estimated to be due to the noise regressors
% is subtracted off.  This may be useful in situations where one wishes to
% treat the denoising as a pre-processing step prior to other analyses 
% of the time-series data.  Further customization of the contents of
% <denoiseddata> is controlled by opt.denoisespec.  If you do not need
% <denoiseddata>, do not assign the <denoiseddata> output when you call 
% GLMdenoisedata.m (this allows us to save on memory usage).
%
% Description of the denoising procedure:
% 1. Determine HRF.  If <hrfmodel> is 'assume', we just use the HRF
%    specified by the user.  If <hrfmodel> is 'optimize', we perform
%    a full fit of the GLM model to the data, optimizing the shape of
%    the HRF.  If <hrfmodel> is 'fir', we do nothing (since full 
%    flexibility in the HRF is allowed for each voxel and each condition).
% 2. Determine cross-validated R^2 values.  We fix the HRF to what is
%    obtained in step 1 and estimate the rest of the GLM model.  Leave-one-
%    run-out cross-validation is performed, and we obtain an estimate of the
%    amount of variance (R^2) that can be predicted by the deterministic 
%    portion of the GLM model (the HRF and the amplitudes).
% 3. Determine noise pool.  This is done by calculating a mean volume (the 
%    mean across all volumes) and then determining the voxels that
%    satisfy the following criteria:
%    (1) The voxels must have sufficient MR signal, that is, the signal
%        intensity in the mean volume must be above a certain threshold
%        (see opt.brainthresh).
%    (2) The voxels must have cross-validated R^2 values that are 
%        below a certain threshold (see opt.brainR2).
%    (3) The voxels must not be listed in opt.brainexclude (which is an
%        optional input that can be specified by the user).
% 4. Determine noise regressors.  This is done by extracting the 
%    time-series data for the voxels in the noise pool, projecting out the
%    polynomial nuisance functions from each time-series, normalizing each
%    time-series to be unit length, and then performing PCA.  The top N
%    PCs from each run (where N is equal to opt.numpcstotry) are selected
%    as noise regressors.  Each regressor is scaled to have a standard
%    deviation of 1; this makes it easier to interpret the weights estimated
%    for the regressors.
% 5. Evaluate different numbers of PCs using cross-validation.  We refit
%    the GLM model to the data (keeping the HRF fixed), systematically varying 
%    the number of PCs from 1 to N.  For each number of PCs, leave-one-run-out 
%    cross-validation is performed.  (Recall that only the deterministic
%    portion of the model is cross-validated; thus, any changes in R^2
%    directly reflect changes in the quality of the amplitude estimates.)
% 6. Choose optimal number of PCs.  To choose the optimal number of PCs,
%    we select a subset of voxels (namely, any voxel that has a cross-validated
%    R^2 value greater than opt.pcR2cutoff (default: 0%) in any of the cases
%    being considered) and then compute the median cross-validated R^2 for these 
%    voxels for different numbers of PCs.  Starting from 0 PCs, we select the 
%    number of PCs that achieves a cross-validation accuracy within opt.pcstop of 
%    the maximum.  (The default for opt.pcstop is 1.05, which means that the
%    chosen number of PCs will be within 5% of the maximum.)
% 7. Fit the final model.  We fit the final GLM model (with the HRF fixed to 
%    that obtained in step 1 and with the number of PCs selected in step 6) 
%    to the data.  Bootstrapping is used to estimate error bars on 
%    amplitude estimates.  (We also fit the final model with no PCs so that
%    we can compare the SNR before and after denoising.)
% 8. Return the de-noised data.  We calculate the component of the data that 
%    is due to the noise regressors and return the original time-series 
%    data with this component subtracted off.  Note that the other components of
%    the model (the hemodynamic responses evoked by the experiment, 
%    the polynomial drift, any extra regressors provided by the user, 
%    model residuals) remain in the de-noised data.  To change this behavior, 
%    please see the input opt.denoisespec.
%
% Figures:
% - "CheckData/CheckMeanStd_runN.png" shows the mean and standard deviation across voxels
%   of each volume in run N.  This allows one to quickly check the sanity of the data, 
%   e.g., with regards to whether weird transient artifacts exist, whether initial 
%   magnetization effects are present in the data (the first few volumes should have
%   been dropped to avoid these effects), whether there are non-finite values in the 
%   data (there should not be any), and with regards to the units of the data 
%   (the data should consist primarily of positive values and in particular, should
%   not be mean-subtracted).
% - "CheckData/CheckDVARS_runN.png" shows the DVARS metric (Power 2012).  Specifically,
%   this is a time-series of the root-mean-square of the difference between pairs of
%   successive volumes.
% - "HRF.png" shows the initial assumed HRF (provided by <hrfknobs>) and the
%   final estimated HRF (as calculated in step 1).  If <hrfmodel> is 'assume',
%   the two plots will be identical.  If <hrfmodel> is 'fir', this figure
%   is not written.
% - "HRFfitmask.png" shows (in white) the mask restricting the voxels that
%   can be chosen for fitting the global HRF.  This figure is written only
%   if <hrfmodel> is 'optimize' and is not written if opt.hrffitmask is 1.
% - "HRFfitvoxels.png" shows (in white) the voxels used to fit the global HRF.
%   This figure is written only if <hrfmodel> is 'optimize' and only if
%   the fitted global HRF is actually chosen for use (in some cases, the
%   initial HRF estimate is chosen; see GLMestimatemodel.m).
% - "PCselection.png" shows for different numbers of PCs, the median 
%   cross-validated R^2 across a subset of voxels (namely, those voxels that 
%   have greater than opt.pcR2cutoff (default: 0%) R^2 for at least one of 
%   the models considered).  The selected number of PCs is circled and 
%   indicated in the title of the figure.
% - "PCscatterN.png" shows a scatter plot of cross-validated R^2 values obtained
%   with no PCs against values obtained with N PCs.  The range of the plot
%   is set to the full range of all R^2 values (across all numbers of PCs).
%   Two different sets of points are plotted.  The first set is shown in green,
%   and this is a set of up to 20,000 voxels randomly selected from the
%   entire pool of voxels.  The second set is shown in red, and this is a set of
%   up to 20,000 voxels randomly selected from the set of voxels that
%   were used to select the number of PC regressors.
% - "MeanVolume.png" shows the mean volume (mean of all volumes).
% - "NoisePool.png" shows (in white) the voxels selected for the noise pool.
% - "NoiseExclude.png" shows (in white) voxels that were excluded by the user
%   from entering the noise pool.  This figure is not written if opt.brainexclude is 0.
% - "PCcrossvalidationN.png" shows cross-validated R^2 values obtained with N PCs.
%   The range is 0% to 100%, and the colormap is nonlinearly scaled to enhance
%   visibility.
% - "PCcrossvalidationscaledN.png" shows cross-validated R^2 values obtained with N PCs.
%   This is just like the previous figures except that the minimum and maximum of the
%   color range is set to the 1st and 99th percentile of R^2 values observed across all 
%   PCs.  The point of this is to make it easier to visualize datasets where R^2 values
%   are low.  The disadvantage is that these figures cannot be directly compared across 
%   different datasets (since the color range may differ).
% - "PCmask.png" shows (in white) the mask restricting the voxels that
%   can be selected for determining the optimal number of PCs.  This figure is 
%   not written if opt.pcR2cutoffmask is 1.
% - "PCvoxels.png" shows (in white) the voxels used to determine the optimal
%   number of PCs.
% - "FinalModel.png" shows R^2 values for the final model (as estimated in
%   step 7).  Note that these R^2 values are not cross-validated.
% - "FinalModel_runN.png" shows R^2 values for the final model separated by 
%   runs.  For example, FinalModel_run01.png indicates R^2 values calculated
%   over the data in run 1.  This might be useful for deciding post-hoc to
%   exclude certain runs from the analysis.
% - "SNRsignal.png" shows the maximum absolute amplitude obtained (the signal).
%   The range is 0 to the 99th percentile of the values.
% - "SNRnoise.png" shows the average amplitude error (the noise).
%   The range is 0 to the 99th percentile of the values.
% - "SNR.png" shows the signal-to-noise ratio.  The range is 0 to 10.
% - "SNRcomparebeforeandafter.png" shows a scatter plot of SNR values obtained
%   with no denoising (0 PCs) against SNR values obtained with denoising (with
%   the selected number of PCs).  The range of the plot is set to the full range
%   of SNR values.  The green and red coloring of the dots is the same as in the
%   "PCscatterN.png" figures.
% - "SNRamountofdatagained.png" is exactly the same as the previous figure except
%   that the y-axis is now converted into the equivalent amount of data gained.  
%   For example, if SNR is boosted from 4 to 8, this is equivalent to having 
%   obtained 300% more data than was actually obtained.
% - "PCmap/PCmap_runN_numO.png" shows the estimated weights for the Oth PC
%   for run N.  The range is -A to A where A is the 99th percentile of the
%   absolute value of all weights across all runs.  The colormap proceeds
%   from blue (negative) to black (0) to red (positive).
%
% Additional information:
% - For additional details on model estimation and quantification of model
% accuracy (R^2), please see the documentation provided in GLMestimatemodel.m.
% - To ensure that <SNRbefore> and <SNRafter> are directly comparable, we first
% average the estimated signal across the no-denoise and denoise cases
% and then divide by the estimated noise from each case.  This has the 
% consequence that <SNRafter> will not be exactly identical to <SNR>.  
% - With regards to the onset-time specification of the experimental design,
% this is implemented by laying down the HRF at the appropriate onset times
% and then using cubic interpolation to obtain the predicted HRF values at the
% times at which data are actually sampled.
%
% History:
% - 2017/02/05: LTD - Adding plots of PC timecourse below pcweight images. and 
%               adding an opt.reconmask. This allows image reconstruction from data that has been vectorized using a mask
%               ie. has no zeros whatsoever. If the logical, 3d mask is load into opt.reconmask
%               and the data is provided as a vector, this will reconstruct all the pretty figures.
% - 2016/09/02: to avoid weird .mat file saving issues, do not save inputs.opt.drawfunction
% - 2016/04/15: add opt.drawfunction
% - 2014/08/01: add opt.wantparametric input (which enables parametric GLM fits).
%               add opt.wantsanityfigures input (which allows user to turn off
%               the sanity-check figures).  add new DVARS sanity-check figures.
%               change the form of the polynomials used in the GLM: the
%               polynomials are now constructed to be orthogonal and unit-length.
%               note that this changes previous behavior and due to numerical
%               issues, the outputs given by this function may be slightly
%               different compared to previous results.
% - 2014/07/14: add <noisepooldirect> input; fix some file- and directory-related
%               issues, making things more platform independent
% - 2013/12/11: rename "global noise regressors" to "noise regressors" in the
%               documentation; perform a quick error-check for non-finite values
%               in the first run of data; omit saving of some of the figures if 
%               opt.numboots is 0; add some new sanity-check figures; fix minor bug
% - 2013/11/18: update documentation; add opt.pccontrolmode; add opt.hrfthresh;
%               use fullfile for compatibility with different platforms;
%               if pcvoxels is empty, we now fallback to using the top 100 voxels;
%               change the colormap scaling for PCcrossvalidationscaled*.png;
%               other various small tweaks;
% - 2013/05/12: allow <design> to specify onset times
% - 2013/05/12: add opt.brainexclude and associated figure
% - 2013/05/12: add SNRbefore and SNRafter fields and associated figures
% - 2013/05/12: add PCcrossvalidationscaledN.png figures
% - 2013/05/12: update to indicate fractional values in design matrix are allowed.
% - 2013/05/12: regressors that are all zero now receive a 0 weight (instead of crashing)
% - 2013/05/12: add pcR2final as an output
% - 2012/12/06: automatically convert data to single format
% - 2012/12/05: fix minor bug (bad HRF estimate caused figure crash)
% - 2012/12/03: *** Tag: Version 1.02 ***. Use faster OLS computation (less
%   error-checking; program execution will halt if design matrix is singular);
%   documentation tweak; minor bug fix.
% - 2012/11/24:
%   - INPUTS: add opt.hrffitmask; opt.pcR2cutoff; opt.pcR2cutoffmask; opt.pcstop; opt.denoisespec; opt.wantpercentbold;
%   - OUTPUTS: add hrffitvoxels, pcvoxels, pcweights, inputs
%   - FIGURES: HRFfitmask.png, HRFfitvoxels.png, PCmask.png, PCvoxels.png, Signal.png, Noise.png, SNR.png, PCmap*.png
%   - hrfmodel can now be 'fir'!
%   - change default of opt.numpcstotry to 20
%   - PC regressors are now scaled to have standard deviation of 1
%   - xval scatter plots now are divided into red and green dots (and up to 20,000 each)
%   - pcselection figure uses a line now (not bar) and a selection circle is drawn
%   - no more skipping of the denoiseddata computation
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
  if isequal(hrfmodel,'fir')
    hrfknobs = 20;
  else
    hrfknobs = normalizemax(getcanonicalhrf(stimdur,tr)');
  end
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
for p=1:length(data)
  if ~isa(data{p},'single')
    fprintf('*** GLMdenoisedata: converting data in run %d to single format (consider doing this before the function call to reduce memory usage). ***\n',p);
    data{p} = single(data{p});
  end
end

% do some error checking
if any(flatten(~isfinite(data{1})))
  fprintf('*** GLMdenoisedata: WARNING: we checked the first run and found some non-finite values (e.g. NaN, Inf). unexpected results may occur due to non-finite values. please fix and re-run GLMdenoisedata. ***\n');
end

% calc
numruns = length(design);
dataclass = class(data{1});  % will always be 'single'
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
numvoxels = prod(xyzsize);

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
if ~isfield(opt,'hrffitmask') || isempty(opt.hrffitmask)
  opt.hrffitmask = 1;
end
if ~isfield(opt,'brainthresh') || isempty(opt.brainthresh)
  opt.brainthresh = [99 0.5];
end
if ~isfield(opt,'brainR2') || isempty(opt.brainR2)
  opt.brainR2 = 0;
end
if ~isfield(opt,'brainexclude') || isempty(opt.brainexclude)
  opt.brainexclude = 0;
end
if ~isfield(opt,'numpcstotry') || isempty(opt.numpcstotry)
  opt.numpcstotry = 20;
end
if ~isfield(opt,'pcR2cutoff') || isempty(opt.pcR2cutoff)
  opt.pcR2cutoff = 0;
end
if ~isfield(opt,'pcR2cutoffmask') || isempty(opt.pcR2cutoffmask)
  opt.pcR2cutoffmask = 1;
end
if ~isfield(opt,'pcstop') || isempty(opt.pcstop)
  opt.pcstop = 1.05;
end
if ~isfield(opt,'pccontrolmode') || isempty(opt.pccontrolmode)
  opt.pccontrolmode = 0;
end
if ~isfield(opt,'numboots') || isempty(opt.numboots)
  opt.numboots = 100;
end
if ~isfield(opt,'denoisespec') || (isempty(opt.denoisespec) && ~iscell(opt.denoisespec))
  opt.denoisespec = '11101';
end
if ~isfield(opt,'wantpercentbold') || isempty(opt.wantpercentbold)
  opt.wantpercentbold = 1;
end
if ~isfield(opt,'hrfthresh') || isempty(opt.hrfthresh)
  opt.hrfthresh = 50;
end
if ~isfield(opt,'noisepooldirect') || isempty(opt.noisepooldirect)
  opt.noisepooldirect = [];
end
if ~isequal(hrfmodel,'fir')
  hrfknobs = normalizemax(hrfknobs);
end
if ~isfield(opt,'wantparametric') || isempty(opt.wantparametric)
  opt.wantparametric = 0;
end
if ~isfield(opt,'wantsanityfigures') || isempty(opt.wantsanityfigures)
  opt.wantsanityfigures = 1;
end
if ~isfield(opt,'drawfunction') || isempty(opt.drawfunction)
  opt.drawfunction = @(vals) vals;
end
if length(opt.maxpolydeg) == 1
  opt.maxpolydeg = repmat(opt.maxpolydeg,[1 numruns]);
end
if ~iscell(opt.extraregressors)
  opt.extraregressors = {opt.extraregressors};
end
if ~iscell(opt.denoisespec)
  opt.denoisespec = {opt.denoisespec};
end

% delete and/or make figuredir
if ~isempty(figuredir)
  if exist(figuredir,'dir')
    assert(rmdir(figuredir,'s'));
  end
  assert(mkdir(figuredir));
  assert(mkdir(fullfile(figuredir,'PCmap')));
  if opt.wantsanityfigures
    assert(mkdir(fullfile(figuredir,'CheckData')));
  end
  figuredir = absolutepath(figuredir);
end

% whether to bypass a lot of the usual GLMdenoise routine
% (since the noise pool and number of PCs are already supplied by the user)
wantbypass = ~isempty(opt.noisepooldirect);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRE-FLIGHT SANITY CHECK

if ~isempty(figuredir) && opt.wantsanityfigures
  fprintf('*** GLMdenoisedata: generating sanity-check figures. ***\n');
  
  % make figure showing mean and std dev of each volume over time
  for p=1:length(data)
  
    % calc
    meants = mean(squish(data{p},dimdata),1);
    stdts = std(squish(data{p},dimdata),[],1);
    dvarts = sqrt(mean(diff(squish(data{p},dimdata),1,2).^2,1));  % RMS of frame-to-frame change

    % make mean+std figure
    figureprep([100 100 1000 300]); hold on;
    errorbar2(1:length(meants),meants,stdts,'v','r-');
    plot(1:length(meants),meants,'r-');
    ax = axis; axis([1-10 length(meants)+10 ax(3:4)]);
    xlabel('TR');
    ylabel('Signal (raw units)');
    title(sprintf('Run %d (mean + std dev of each volume)',p));
    figurewrite(sprintf('CheckMeanStd_run%02d',p),[],[],fullfile(figuredir,'CheckData'));

    % make dvars figure
    figureprep([100 100 1000 300]); hold on;
    plot(dvarts,'r-');
    ax = axis; axis([1-10 length(dvarts)+10 ax(3:4)]);
    xlabel('Volume number');
    ylabel('RMS of difference image (raw units)');
    title(sprintf('Run %d (DVARS)',p));
    figurewrite(sprintf('CheckDVARS_run%02d',p),[],[],fullfile(figuredir,'CheckData'));

  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DETERMINE HRF

% if 'optimize', perform full-fit to determine HRF
switch hrfmodel
case 'optimize'
  fprintf('*** GLMdenoisedata: performing full fit to estimate global HRF. ***\n');
  fullfit = GLMestimatemodel(design,data,stimdur,tr,hrfmodel,hrfknobs,0,opt,[],2);
  hrf = fullfit.modelmd{1};
  hrffitvoxels = fullfit.hrffitvoxels;
  clear fullfit;

% if 'assume', the HRF is provided by the user
case 'assume'
  hrf = hrfknobs;
  hrffitvoxels = [];

% if 'fir', do nothing
case 'fir'
  hrffitvoxels = [];
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALCULATE CROSS-VALIDATED R^2 VALUES

if ~wantbypass

  % perform cross-validation to determine R^2 values
  fprintf('*** GLMdenoisedata: performing cross-validation to determine R^2 values. ***\n');
  switch hrfmodel
  case {'optimize' 'assume'}
    xvalfit = GLMestimatemodel(design,data,stimdur,tr,'assume',hrf,-1,opt,[],1);
  case 'fir'
    xvalfit = GLMestimatemodel(design,data,stimdur,tr,'fir',hrfknobs,-1,opt,[],1);
  end
  pcR2 = xvalfit.R2;
  clear xvalfit;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DETERMINE NOISE POOL AND CALCULATE NOISE REGRESSORS

% calculate mean volume
volcnt = cellfun(@(x) size(x,dimtime),data);
meanvol = reshape(catcell(2,cellfun(@(x) squish(mean(x,dimtime),dimdata),data,'UniformOutput',0)) ...
                  * (volcnt' / sum(volcnt)),[xyzsize 1]);

if ~wantbypass

  % determine noise pool
  fprintf('*** GLMdenoisedata: determining noise pool. ***\n');
  thresh = prctile(meanvol(:),opt.brainthresh(1))*opt.brainthresh(2);  % threshold for non-brain voxels 
  bright = meanvol > thresh;                                           % logical indicating voxels that are bright (brain voxels)
  badxval = pcR2 < opt.brainR2;                                        % logical indicating voxels with poor cross-validation accuracy
  noisepool = bright & badxval & ~opt.brainexclude;                    % logical indicating voxels that satisfy all criteria

else

  % the noise pool is specified by the user
  noisepool = opt.noisepooldirect{1};

end
  
% determine noise regressors
fprintf('*** GLMdenoisedata: calculating noise regressors. ***\n');
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
  u = bsxfun(@rdivide,u,std(u,[],1));  % scale so that std is 1
  pcregressors{p} = cast(u,dataclass);

end
clear temp len u s v;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PERTURB THE NOISE REGRESSORS (IF REQUESTED)

switch opt.pccontrolmode

% do nothing (this is the default)
case 0

% phase-scramble each regressor
case 1

  % for each run
  for p=1:length(pcregressors)
  
    % for each regressor
    for q=1:size(pcregressors{p},2)
    
      % the original regressor
      temp = pcregressors{p}(:,q);

      % a sample of white noise
      temp2 = randn(size(temp));

      % new regressor has the amplitude spectrum of the original regressor,
      % but the phase spectrum of the white noise
      pcregressors{p}(:,q) = real(ifft(abs(fft(temp,[],1)) .* exp(j * angle(fft(temp2,[],1)))));

    end

  end
  clear temp temp2;

% shuffle regressors across runs (ensuring none match up to the original assignment)
case 2

  % repeat until we have a satisfactory assignment
  while 1
  
    % generate a random permutation
    temp = permutedim(1:length(pcregressors));
    
    % if none matched the original assignment, we are done
    if ~any(temp == (1:length(pcregressors)))
      break;
    end

  end

  % shuffle the regressors
  pcregressors = pcregressors(temp);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADD NOISE REGRESSORS INTO MODEL AND CHOOSE OPTIMAL NUMBER

if ~wantbypass

  % perform cross-validation with increasing number of noise regressors
  for p=1:opt.numpcstotry
    fprintf('*** GLMdenoisedata: performing cross-validation with %d PCs. ***\n',p);
    opt2 = opt;
    for q=1:numruns
      opt2.extraregressors{q} = cat(2,opt2.extraregressors{q},pcregressors{q}(:,1:p));
    end
    opt2.wantpercentbold = 0;  % no need to do this, so optimize for speed
    switch hrfmodel
    case {'optimize' 'assume'}
      xvalfit = GLMestimatemodel(design,data,stimdur,tr,'assume',hrf,-1,opt2,[],1);
    case 'fir'
      xvalfit = GLMestimatemodel(design,data,stimdur,tr,'fir',hrfknobs,-1,opt2,[],1);
    end
    pcR2 = cat(dimdata+1,pcR2,xvalfit.R2);
  end
  clear xvalfit;

  % prepare to select optimal number of PCs
  temp = squish(pcR2,dimdata);  % voxels x 1+pcs
  pcvoxels = any(temp > opt.pcR2cutoff,2) & squish(opt.pcR2cutoffmask,dimdata);  % if pcR2cutoffmask is 1, this still works
  if ~any(pcvoxels)
    warning(['No voxels passed the threshold for the selection of the number of PCs. ' ...
             'We fallback to simply using the top 100 voxels (i.e. compute each voxel''s maximum ' ...
             'cross-validation accuracy under any number of PCs and then choose the top 100 voxels.']);
    if isequal(opt.pcR2cutoffmask,1)
      ix = 1:size(temp,1);
    else
      ix = find(squish(opt.pcR2cutoffmask,dimdata));
    end
    pcvoxels = logical(zeros(size(temp,1),1));
    temp2 = max(temp(ix,:),[],2);  % max cross-validation for each voxel (within mask)
    [d,ix2] = sort(temp2,'descend');
    pcvoxels(ix(ix2(1:min(100,length(ix2))))) = 1;
  end
  xvaltrend = median(temp(pcvoxels,:),1);

  % choose number of PCs
  chosen = 0;  % this is the fall-back
  if opt.pcstop <= 0
    chosen = -opt.pcstop;  % in this case, the user decides
  else
    curve = xvaltrend - xvaltrend(1);  % this is the performance curve that starts at 0 (corresponding to 0 PCs)
    mx = max(curve);                   % store the maximum of the curve
    best = -Inf;                       % initialize (this will hold the best performance observed thus far)
    for p=0:opt.numpcstotry
  
      % if better than best so far
      if curve(1+p) > best
    
        % record this number of PCs as the best
        chosen = p;
        best = curve(1+p);
      
        % if we are within opt.pcstop of the max, then we stop.
        if best*opt.pcstop >= mx
          break;
        end
      
      end
  
    end
  end

  % record the number of PCs
  pcnum = chosen;
  fprintf('*** GLMdenoisedata: selected number of PCs is %d. ***\n',pcnum);

else

  % the number of PCs is specified by the user
  pcnum = opt.noisepooldirect{2};
  fprintf('*** GLMdenoisedata: user-specified number of PCs is %d. ***\n',pcnum);
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIT FINAL MODEL AND PREPARE OUTPUT

% fit final model (NO DENOISING)
opt2 = opt;
opt2.wantpercentbold = 0;  % do not do the conversion yet.  we will do it ourselves below.
fprintf('*** GLMdenoisedata: fitting final model (no denoising, for comparison purposes). ***\n');
switch hrfmodel
case {'optimize' 'assume'}
  resultsTEMP = GLMestimatemodel(design,data,stimdur,tr,'assume',hrf,opt.numboots,opt2);
case 'fir'
  resultsTEMP = GLMestimatemodel(design,data,stimdur,tr,'fir',hrfknobs,opt.numboots,opt2);
end

% save some results
signal_nodenoise = resultsTEMP.signal;
noise_nodenoise = resultsTEMP.noise;
clear resultsTEMP;

% fit final model (WITH DENOISING)
opt2 = opt;
for q=1:numruns
  opt2.extraregressors{q} = cat(2,opt2.extraregressors{q},pcregressors{q}(:,1:pcnum));
end
opt2.wantpercentbold = 0;  % do not do the conversion yet.  we will do it ourselves below.
fprintf('*** GLMdenoisedata: fitting final model (with denoising). ***\n');
switch hrfmodel     % note that we will use the rawdesign field from the cache output for the parametric fits
case {'optimize' 'assume'}
  [results,cache] = GLMestimatemodel(design,data,stimdur,tr,'assume',hrf,opt.numboots,opt2);
case 'fir'
  [results,cache] = GLMestimatemodel(design,data,stimdur,tr,'fir',hrfknobs,opt.numboots,opt2);
end

% prepare additional outputs
results.hrffitvoxels = hrffitvoxels;  % note that this overrides the existing entry in results
results.meanvol = meanvol;
results.noisepool = noisepool;
results.pcregressors = pcregressors;
if ~wantbypass
  results.pcR2 = pcR2;
  results.pcvoxels = reshape(pcvoxels,[xyzsize 1]);
end
results.pcnum = pcnum;
clear meanvol noisepool pcregressors pcR2 pcnum;

% prepare SNR comparison
amp = (results.signal + signal_nodenoise)/2;
results.SNRbefore = amp ./ noise_nodenoise;
results.SNRafter = amp ./ results.noise;
clear amp signal_nodenoise noise_nodenoise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETRIC FITS AND ERROR ESTIMATES

if opt.wantparametric

  fprintf('*** GLMdenoisedata: calculating parametric fits and error estimates. ***\n');

  % construct design matrix
  temp = {};
  for p=1:numruns
    numtime = size(data{p},dimtime);
    temp{p} = cat(2,constructpolynomialmatrix(numtime,0:opt.maxpolydeg(p)), ...
                    opt.extraregressors{p}, ...
                    results.pcregressors{p}(:,1:results.pcnum));
  end
  results.parametric.designmatrix = [catcell(1,cache.rawdesign) blkdiag(temp{:})];  % time x regressors

  % estimate parameters using OLS
  results.parametric.parameters = ...
    mtimescell(olsmatrix2(results.parametric.designmatrix), ...
               cellfun(@(x) squish(x,dimdata)',data,'UniformOutput',0));  % regressors x voxels

  % compute sum of the squares of the residuals (e.g. sum(resid.^2))
  sumsq = sum((results.parametric.designmatrix * results.parametric.parameters - ...
               catcell(1,cellfun(@(x) squish(x,dimdata)',data,'UniformOutput',0))).^2,1);  % 1 x voxels

  % estimate noise variance (e.g. sum(resid.^2) / (n - rank(X)))
  results.parametric.noisevar = ...  % 1 x voxels
    sumsq ./ (size(results.parametric.designmatrix,1) - rank(results.parametric.designmatrix));

  % calculate standard error on parameters (e.g. sqrt(sigma^2 * inv(X'*X)))
  good = ~all(results.parametric.designmatrix==0,1);
  X = results.parametric.designmatrix(:,good);
  results.parametric.parametersse = ...  % regressors x voxels
    sqrt(bsxfun(@times,results.parametric.noisevar, ...
                copymatrix(zeros(size(results.parametric.designmatrix,2),1),good,diag(inv(X'*X)))));
  
  % clean up
  clear temp sumsq good X;
  
  % prepare for output
  results.parametric.parameters =   reshape(results.parametric.parameters',  [xyzsize size(results.parametric.parameters,1)]);
  results.parametric.noisevar =     reshape(results.parametric.noisevar,     [xyzsize 1]);
  results.parametric.parametersse = reshape(results.parametric.parametersse',[xyzsize size(results.parametric.parametersse,1)]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALCULATE DENOISED DATA AND PCWEIGHTS

fprintf('*** GLMdenoisedata: calculating denoised data and PC weights. ***\n');

% for each run, perform regression to figure out the various contributions
denoiseddata = cell(length(opt.denoisespec),numruns);
for p=1:numel(denoiseddata)
  denoiseddata{p} = cast(denoiseddata{p},dataclass);
end
results.pcweights = zeros([numvoxels results.pcnum numruns],dataclass);
for p=1:numruns

  % calc
  numtime = size(data{p},dimtime);

  % calculate signal contribution
  modelcomponent = GLMpredictresponses(results.modelmd,{design{p}},tr,numtime,dimdata);
  modelcomponent = modelcomponent{1};  % X x Y x Z x T

  % prepare polynomial regressors
  polymatrix = constructpolynomialmatrix(numtime,0:opt.maxpolydeg(p));
  numpoly = size(polymatrix,2);

  % prepare other regressors
  exmatrix = opt.extraregressors{p};
  numex = size(exmatrix,2);

  % prepare noise regressors
  pcmatrix = results.pcregressors{p}(:,1:results.pcnum);
  numpc = size(pcmatrix,2);

  % estimate weights
  h = olsmatrix2(cat(2,polymatrix,exmatrix,pcmatrix))*squish(data{p} - modelcomponent,dimdata)';  % parameters x voxels

  % record weights on noise regressors
  results.pcweights(:,:,p) = h(numpoly+numex+(1:numpc),:)';
  
  % construct time-series
  polycomponent = reshape((polymatrix*h(1:numpoly,:))',[xyzsize numtime]);
  if numex == 0
    excomponent = zeros([xyzsize numtime],dataclass);
  else
    excomponent = reshape((exmatrix*h(numpoly+(1:numex),:))',[xyzsize numtime]);
  end
  if numpc == 0
    pccomponent = zeros([xyzsize numtime],dataclass);
  else
    pccomponent = reshape((pcmatrix*h(numpoly+numex+(1:numpc),:))',[xyzsize numtime]);
  end
  residcomponent = data{p} - (modelcomponent + polycomponent + excomponent + pccomponent);
  
  % construct denoised data (but don't bother if the user doesn't want it)
  if nargout > 1
    for q=1:length(opt.denoisespec)
      denoiseddata{q,p} = bitget(bin2dec(opt.denoisespec{q}),5) * modelcomponent + ...
                         bitget(bin2dec(opt.denoisespec{q}),4) * polycomponent + ...
                         bitget(bin2dec(opt.denoisespec{q}),3) * excomponent + ...
                         bitget(bin2dec(opt.denoisespec{q}),2) * pccomponent + ...
                         bitget(bin2dec(opt.denoisespec{q}),1) * residcomponent;
    end
  end

end

% clean up
clear modelcomponent h polycomponent excomponent pccomponent residcomponent;

% prepare for output
results.pcweights = reshape(results.pcweights,[xyzsize results.pcnum numruns]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE ADDITIONAL OUTPUTS

% return all the inputs (except for the data) in the output.
% also, include a new field 'datasize'.
results.inputs.design = design;
results.inputs.datasize = cellfun(@(x) size(x),data,'UniformOutput',0);
results.inputs.stimdur = stimdur;
results.inputs.tr = tr;
results.inputs.hrfmodel = hrfmodel;
results.inputs.hrfknobs = hrfknobs;
results.inputs.opt = opt;
results.inputs.opt = rmfield(results.inputs.opt,'drawfunction');  % anonymous functions sometimes cause major file-saving problems, so just omit
results.inputs.figuredir = figuredir;

if ~wantbypass
  % prepare pcR2final
  results.pcR2final = subscript(results.pcR2,[repmat({':'},[1 dimdata]) {1+results.pcnum}]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONVERT TO % BOLD CHANGE

if opt.wantpercentbold
  fprintf('*** GLMdenoisedata: converting to percent BOLD change. ***\n');
  con = 1./abs(results.meanvol) * 100;
  switch hrfmodel
  case 'fir'
    for p=1:size(results.models,dimdata+3)  % ugly to save MEMORY
      if dimdata==3
        results.models(:,:,:,:,:,p) = bsxfun(@times,results.models(:,:,:,:,:,p),con);
      else
        results.models(:,:,:,p) = bsxfun(@times,results.models(:,:,:,p),con);
      end
    end
    results.modelmd = bsxfun(@times,results.modelmd,con);
    results.modelse = bsxfun(@times,results.modelse,con);
  case {'assume' 'optimize'}
    for p=1:size(results.models{2},dimdata+2)  % ugly to save MEMORY
      if dimdata==3
        results.models{2}(:,:,:,:,p) = bsxfun(@times,results.models{2}(:,:,:,:,p),con);
      else
        results.models{2}(:,:,p) = bsxfun(@times,results.models{2}(:,:,p),con);
      end
    end
    results.modelmd{2} = bsxfun(@times,results.modelmd{2},con);
    results.modelse{2} = bsxfun(@times,results.modelse{2},con);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERATE FIGURES

if ~isempty(figuredir)
  fprintf('*** GLMdenoisedata: generating figures. ***\n');

  % make figure showing HRF
  if ~isequal(hrfmodel,'fir') && length(hrfknobs) > 1
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
  end
  
   %create the correctly sized mask for reconstructing 3d images from vecs
    if isfield(opt, 'reconmask')
        mask_vec = reshape(opt.reconmask,[],1);
        x_dim = size(opt.reconmask,1);
        y_dim = size(opt.reconmask,2);
        z_dim = size(opt.reconmask,3);
    end
    
  
  % write out image showing HRF fit voxels
  if isequal(hrfmodel,'optimize') && ~isempty(results.hrffitvoxels)
        if isfield(opt, 'reconmask')
            mask_vec = reshape(opt.reconmask,[],1);
            hrffit_recon = mask_vec;
            hrffit_recon(mask_vec ~=0) = results.hrffitvoxels;
            hrffit_recon = single(reshape(hrffit_recon,x_dim,y_dim,z_dim));
            results.hrffitvoxels = hrffit_recon;
        end
    imwrite(uint8(255*makeimagestack(opt.drawfunction(results.hrffitvoxels),[0 1])),gray(256),fullfile(figuredir,'HRFfitvoxels.png'));
  end

  if ~wantbypass

    % make figure illustrating selection of number of PCs
    figureprep([100 100 400 400]); hold on;
    plot(0:opt.numpcstotry,xvaltrend,'r.-');
    set(scatter(results.pcnum,xvaltrend(1+results.pcnum),100,'ro'),'LineWidth',2);
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
      scattersparse(temp(:,1),temp(:,1+p),20000,0,36,'g','.');
      scattersparse(temp(pcvoxels,1),temp(pcvoxels,1+p),20000,0,36,'r','.');
      axis([rng rng]); axissquarify; axis([rng rng]); 
      straightline(0,'h','y-');
      straightline(0,'v','y-');
      xlabel('Cross-validated R^2 (0 PCs)');
      ylabel(sprintf('Cross-validated R^2 (%d PCs)',p));
      title(sprintf('Number of PCs = %d',p));
      figurewrite(sprintf('PCscatter%02d',p),[],[],figuredir);
    end

  end

  % write out image showing mean volume (of first run)
  if isfield(opt, 'reconmask')
        mean_vol_recon = mask_vec;
        mean_vol_recon(mask_vec ~=0) = results.meanvol;
        mean_vol_recon = reshape(mean_vol_recon, x_dim,y_dim,z_dim);
        results.meanvol=mean_vol_recon;
  end
  imwrite(uint8(255*makeimagestack(opt.drawfunction(results.meanvol),1)),gray(256),fullfile(figuredir,'MeanVolume.png'));

  % write out image showing noise pool
  if ~isempty(results.noisepool)  % in certain degenerate cases, this might happen
        if isfield(opt, 'reconmask')
            noisepool_recon = mask_vec;
            noisepool_recon(mask_vec ~=0) = results.noisepool;
            noisepool_recon = reshape(noisepool_recon, x_dim,y_dim,z_dim);
            results.noisepool=noisepool_recon;
        end
    imwrite(uint8(255*makeimagestack(opt.drawfunction(results.noisepool),[0 1])),gray(256),fullfile(figuredir,'NoisePool.png'));
  end
  
  % write out image showing voxels excluded from noise pool
  if ~isequal(opt.brainexclude,0)
        if isfield(opt, 'reconmask')
            temp_recon = mask_vec;
            temp_recon(mask_vec ~=0) = opt.brainexclude;
            temp_recon = reshape(temp_recon, x_dim,y_dim,z_dim);
            opt.brainexclude=temp_recon;
        end
    imwrite(uint8(255*makeimagestack(opt.drawfunction(opt.brainexclude),[0 1])),gray(256),fullfile(figuredir,'NoiseExclude.png'));
  end

  % write out image showing voxel mask for HRF fitting
  if isequal(hrfmodel,'optimize') && ~isequal(opt.hrffitmask,1)
       if isfield(opt, 'reconmask')
            temp_recon = mask_vec;
            temp_recon(mask_vec ~=0) = opt.hrffitmask;
            temp_recon = reshape(temp_recon, x_dim,y_dim,z_dim);
            opt.hrffitmask=temp_recon;
        end
    imwrite(uint8(255*makeimagestack(opt.drawfunction(opt.hrffitmask),[0 1])),gray(256),fullfile(figuredir,'HRFfitmask.png'));
  end

  % define a function that will write out R^2 values to an image file
  imfun = @(results,filename) ...
    imwrite(uint8(255*makeimagestack(opt.drawfunction(signedarraypower(results/100,0.5)),[0 1])),hot(256),filename);

  if ~wantbypass

    % write out image showing voxel mask for PC selection
    if ~isequal(opt.pcR2cutoffmask,1)
            if isfield(opt, 'reconmask')
                temp_recon = mask_vec;
                temp_recon(mask_vec ~=0) = opt.pcR2cutoffmask;
                temp_recon = reshape(temp_recon, x_dim,y_dim,z_dim);
                opt.pcR2cutoffmask=temp_recon;
            end
      imwrite(uint8(255*makeimagestack(opt.drawfunction(opt.pcR2cutoffmask),[0 1])),gray(256),fullfile(figuredir,'PCmask.png'));
    end
  
    % write out image showing the actual voxels used for PC selection
    if isfield(opt, 'reconmask')
            temp_recon = mask_vec;
            temp_recon(mask_vec ~=0) = results.pcvoxels;
            temp_recon = reshape(temp_recon, x_dim,y_dim,z_dim);
            results.pcvoxels=temp_recon;
    end
    imwrite(uint8(255*makeimagestack(opt.drawfunction(results.pcvoxels),[0 1])),gray(256),fullfile(figuredir,'PCvoxels.png'));

    % figure out bounds for the R^2 values
    bounds = prctile(results.pcR2(:),[1 99]);
    if bounds(1)==bounds(2)  % a hack to avoid errors in normalization
      bounds(2) = bounds(1) + 1;
    end

    % define another R^2 image-writing function
    imfunB = @(results,filename) ...
      imwrite(uint8(255*makeimagestack(opt.drawfunction(signedarraypower(normalizerange(results,0,1,bounds(1),bounds(2)),0.5)),[0 1])),hot(256),filename);

    % write out cross-validated R^2 for the various numbers of PCs
    for p=1:size(results.pcR2,dimdata+1)
      temp = subscript(results.pcR2,[repmat({':'},[1 dimdata]) {p}]);
      if isfield(opt, 'reconmask')
                temp_recon = mask_vec;
                temp_recon(mask_vec ~=0) = temp;
                temp_recon = reshape(temp_recon, x_dim,y_dim,z_dim);
                temp=temp_recon;
            end
      feval(imfun,temp,fullfile(figuredir,sprintf('PCcrossvalidation%02d.png',p-1)));
      feval(imfunB,temp,fullfile(figuredir,sprintf('PCcrossvalidationscaled%02d.png',p-1)));
    end

  end

  % write out overall R^2 for final model
  if isfield(opt, 'reconmask')
        temp_recon = mask_vec;
        temp_recon(mask_vec ~=0) = results.R2;
        temp_recon = reshape(temp_recon, x_dim,y_dim,z_dim);
        results.R2=temp_recon;
  end
  feval(imfun,results.R2,fullfile(figuredir,sprintf('FinalModel.png')));

  % write out R^2 separated by runs for final model
  for p=1:size(results.R2run,dimdata+1)
    temp = subscript(results.R2run,[repmat({':'},[1 dimdata]) {p}]);
    if isfield(opt, 'reconmask')
            temp_recon = mask_vec;
            temp_recon(mask_vec ~=0) = temp;
            temp_recon = reshape(temp_recon, x_dim,y_dim,z_dim);
            temp=temp_recon;
        end
    feval(imfun,temp,fullfile(figuredir,sprintf('FinalModel_run%02d.png',p)));
  end
  
  % write out signal, noise, and SNR
  if isfield(opt, 'reconmask')
        temp_recon = mask_vec;
        temp_recon(mask_vec ~=0) = results.signal;
        temp_recon = reshape(temp_recon, x_dim,y_dim,z_dim);
        results.signal=temp_recon;
  end
  imwrite(uint8(255*makeimagestack(opt.drawfunction(results.signal),[0 prctile(results.signal(:),99)])),hot(256),fullfile(figuredir,'SNRsignal.png'));
  if opt.numboots ~= 0
                if isfield(opt, 'reconmask')
                    temp_recon = mask_vec;
                    temp_recon(mask_vec ~=0) = results.noise;
                    temp_recon = reshape(temp_recon, x_dim,y_dim,z_dim);
                    results.noise=temp_recon;
                end
    imwrite(uint8(255*makeimagestack(opt.drawfunction(results.noise),[0 max(eps,prctile(results.noise(:),99))])),hot(256),fullfile(figuredir,'SNRnoise.png'));
                if isfield(opt, 'reconmask')
                    temp_recon = mask_vec;
                    temp_recon(mask_vec ~=0) = results.SNR;
                    temp_recon = reshape(temp_recon, x_dim,y_dim,z_dim);
                    results.SNR=temp_recon;
                end
    imwrite(uint8(255*makeimagestack(opt.drawfunction(results.SNR),[0 10])),hot(256),fullfile(figuredir,'SNR.png'));
  end
  
  % write out SNR comparison figures (first figure)
  if opt.numboots ~= 0
    rng = [min([results.SNRbefore(:); results.SNRafter(:)]) max([results.SNRbefore(:); results.SNRafter(:)])];
    if ~all(isfinite(rng))  % hack to deal with cases of no noise estimate
      rng = [0 1];
    end
    figureprep([100 100 500 500]); hold on;
    scattersparse(results.SNRbefore(:),results.SNRafter(:),20000,0,36,'g','.');
    if ~wantbypass
      scattersparse(results.SNRbefore(pcvoxels),results.SNRafter(pcvoxels),20000,0,36,'r','.');
    end
    axis([rng rng]); axissquarify; axis([rng rng]);
    xlabel('SNR (before denoising)');
    ylabel('SNR (after denoising)');
    figurewrite(sprintf('SNRcomparebeforeandafter'),[],[],figuredir);
  end
  
  % write out SNR comparison figures (second figure)
  if opt.numboots ~= 0
    datagain = ((results.SNRafter./results.SNRbefore).^2 - 1) * 100;
    figureprep([100 100 500 500]); hold on;
    scattersparse(results.SNRbefore(:),datagain(:),20000,0,36,'g','.');
    if ~wantbypass
      scattersparse(results.SNRbefore(pcvoxels),datagain(pcvoxels),20000,0,36,'r','.');
    end
    ax = axis; axis([rng ax(3:4)]);
    xlabel('SNR (before denoising)');
    ylabel('Equivalent amount of data gained (%)');
    figurewrite(sprintf('SNRamountofdatagained'),[],[],figuredir);
  end
  
  % write out maps of pc weights
  thresh = prctile(abs(results.pcweights(:)),99);
  for p=1:size(results.pcweights,dimdata+1)
    for q=1:size(results.pcweights,dimdata+2)
      temp = subscript(results.pcweights,[repmat({':'},[1 dimdata]) {p} {q}]);
                    if isfield(opt, 'reconmask')
                        temp_recon = mask_vec;
                        temp_recon(mask_vec ~=0) = temp;
                        temp_recon = reshape(temp_recon, x_dim,y_dim,z_dim);
                        temp=temp_recon;
                    end
                    
                    %Go ahead and create the map as it was originaly done
                    temp = (255*makeimagestack(opt.drawfunction(temp),[-thresh thresh]));
                    
                    time_plot_dim = size(results.pcregressors{1,q}(:,p),1); 
                    %Make sure the x axis is going to have the right
                    %length
                    
                    
                    %Create a 3x3 plot space, use the majority for the
                    %weightmap, only a the last row for weight plot.
                    
                    subplot(4,3,[1, 2, 3, 4, 5, 6, 7, 8, 9]); imshow(temp,cmapsign(256));
                    subplot(4,3,[10, 11, 12]); plot(results.pcregressors{1,q}(:,p)); xlim([0, time_plot_dim]);
                    
                    %save everything to the right place
                    label = strcat('PCmap_run',num2str(q), '_num' , num2str(p));
                    print(fullfile(figuredir, 'PCmap', label), '-dpng');
    end
  end

end
