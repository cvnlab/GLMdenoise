function f = catcell(dim,m)

% function f = catcell(dim,m)
%
% <dim> is the dimension to concatenate along
% <m> is a cell matrix
%
% simply return cat(dim,m{:}).  this function is useful because 
% MATLAB doesn't provide an easy way to apply "{:}" to an 
% arbitrary matrix.
%
% example:
% isequal(catcell(2,{1 2 3}),[1 2 3])

f = cat(dim,m{:});
