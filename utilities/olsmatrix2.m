function f = olsmatrix2(X)

% function f = olsmatrix2(X)
%
% <X> is samples x parameters
%
% what we want to do is to perform OLS regression using <X>
% and obtain the parameter estimates.  this is accomplished
% by inv(X'*X)*X'*y = f*y where y is the data (samples x cases).
%
% what this function does is to return <f> which has dimensions
% parameters x samples.
%
% if any warning messages are produced by the inversion process, then we die.
% this is a conservative strategy that ensures that the regression is 
% well-behaved (i.e. has a unique, finite solution).
%
% note that no scale normalization of the regressor columns is performed.
% also, note that we use \ to perform the inversion.
%
% see also olsmatrix.m.

lastwarn('');
f = (X'*X)\X';
assert(isempty(lastwarn),lastwarn);
