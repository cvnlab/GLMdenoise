function responses = GLMpredictresponses(model,design,dimdata)

% function responses = GLMpredictresponses(model,design,dimdata)
%
% <model> is one of the following:
%   A where A is X x Y x Z x conditions x time with the timecourse of the 
%     response of each voxel to each condition.  XYZ can be collapsed.
%   {B C} where B is time x 1 with the HRF that is common to all voxels and
%     all conditions and C is X x Y x Z x conditions with the amplitude of the 
%     response of each voxel to each condition
%   Note that in both of these cases, the first time point is assumed to be 
%   coincident with condition onset.
% <design> is the experimental design with dimensions time x conditions.
%   Each column should be zeros except for ones indicating condition onsets.
%   Can be a cell vector whose elements correspond to different runs.
%   Different runs can have different numbers of time points.
% <dimdata> indicates the dimensionality of the voxels.
%   A value of 3 indicates X x Y x Z, and a value of 1 indicates XYZ.
%
% Given the model (<model>) and the experimental design (<design>),
% compute the predicted time-series response.
%
% Return:
% <responses> as X x Y x Z x time or a cell vector of elements that are 
%   each X x Y x Z x time.  The format of <responses> will mirror the format
%   of <design> (for example, if <design> is a cell vector, <responses> will 
%   also be a cell vector).
%
% History:
% - 2012/11/2 - Initial version.

% calc
iscellcase = iscell(design);
dimtime = dimdata + 2;
if iscell(model)
  xyzsize = sizefull(model{2},dimdata);
else
  xyzsize = sizefull(model,dimdata);
end

% make cell
if ~iscell(design)
  design = {design};
end

% loop over runs
responses = {};
for p=1:length(design)

  % number of time points in the run
  ntime = size(design{p},1);

  % case of shared HRF
  if iscell(model)
    
    % convolve with HRF
    temp = conv2(full(design{p}),model{1});  % make full just in case design is sparse

    % extract desired subset of time-series
    temp = temp(1:ntime,:);  % time x conditions

    % weight by the amplitudes
    responses{p} = reshape((temp * squish(model{2},dimdata)')',[xyzsize ntime]);  % X x Y x Z x time
  
  % case of individual timecourses
  else

    % length of each timecourse (L)
    len = size(model,dimtime);
    
    % expand design matrix using delta functions
    temp = constructstimulusmatrices(design{p}',0,len-1,0);  % time x L*conditions

    % weight design matrix by the timecourses
    responses{p} = reshape((temp * squish(permute(squish(model,dimdata),[3 2 1]),2))',[xyzsize ntime]);  % X x Y x Z x time

  end

end

% undo cell if necessary
if ~iscellcase
  responses = responses{1};
end
