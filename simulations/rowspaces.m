function [rowbasis,nullspace,dimsfound] = rowspaces(A,dimensions)
% rowspaces: orthogonal basis vectors for the null space and row space of A
% usage: [rowbasis,nullspace,dimsfound] = rowspaces(A)
% usage: [rowbasis,nullspace,dimsfound] = rowspaces(A,dimensions)
%
% rowspaces is a bit faster than null or orth, when
% you just need raw speed.
%
% arguments: (input)
%  A   - nxp array A
%
%  dimensions - (OPTIONAL) Expected number of dimensions in
%        the data. If the expected dimensionality is unknown
%        in advance, then leave this argument out.
%
% arguments: (output)
%  rowbasis  - an orthogonal basis for the rows of A
%
%  nullspace - an orthogonal basis for the nullspace of the
%        rows of A
%
% 
% Example usage:
%   
% A = rand(2,4)
% A =
%       0.7015       0.1056       0.2184      0.37002
%       0.21091      0.32589      0.86503      0.80694
% 
% [rs,ns] = rowspaces(A)
% rs =
%      -0.1694     -0.26175     -0.69478     -0.64813
%       0.96701    -0.055813     -0.24623     0.033755
%
% ns =
%      0.097927     -0.79885      0.52868     -0.26971
%     -0.16309     -0.53871     -0.42088      0.71137
%
% See also: null, orth
%
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 2/24/07

% Use a pivoted qr to get the bases.
% Really, I don't care about the pivot sequence,
% but the use of pivoting should make it more
% stable.
[Q,R,E] = qr(A'); %#ok

dr = diag(R);

tol = 1e4*eps(max(abs(dr)));
k = abs(dr)>tol;
k = [k;false(size(A,2)-length(k),1)];

% How many dimensions did we find in this data?
dimsfound = sum(k);

% expected dimensions supplied?
if (nargin<2) || isempty(dimensions)
  dimensions = dimsfound;
elseif dimensions>size(A,2)
  error('dimensions is incompatible with the size of A')
end

% Throw a warning only if the expected
% dimensionanlity was supplied.
if dimsfound<dimensions
  warning('ROWSPACES:dimensions', ...
    'Unexpectedly degenerate row dimensionality')
elseif dimsfound>dimensions
  warning('ROWSPACES:dimensions', ...
    'Unexpectedly overabundant row dimensionality')
end

% we still return the expected sizes of the arrays
rowbasis = Q(:,1:dimensions)';
nullspace = Q(:,(dimensions+1):end)';


