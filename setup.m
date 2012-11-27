% This script downloads the example dataset (if it has not already
% been downloaded) and adds GLMdenoise to the MATLAB path.

% Download the dataset if it has not already been downloaded.
% (If you like, you may manually download the dataset from:
%    http://kendrickkay.net/GLMdenoise/exampledataset.mat)
files = {'exampledataset.mat'};
for p=1:length(files)
  if ~exist(files{p},'file')
    fprintf('Downloading %s (please be patient).\n',files{p});
    assert(unix(sprintf('curl -O http://kendrickkay.net/GLMdenoise/%s',files{p})) == 0);
    fprintf('Downloading is done!\n');
  end
end
clear p files;

% Add GLMdenoise to the MATLAB path (in case the user has not already done so).
path0 = strrep(which('setup.m'),'/GLMdenoise/setup.m','');
addpath([path0 '/GLMdenoise']);
addpath([path0 '/GLMdenoise/utilities']);
clear path0;
