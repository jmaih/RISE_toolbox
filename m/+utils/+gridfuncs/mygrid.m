function [Regimes,unique_regimes]=mygrid(v)
% mygrid creates a grid of points
%
% Syntax
% -------
% ::
%
%   [Regimes,unique_regimes]=mygrid(v)
%
% Inputs
% -------
%
% - **v** [scalar|vector] : n x 1 or 1 x n each entry in the vector
%   represents the number of elements along that dimension
%
% Outputs
% --------
%
% - **Regimes** [matrix|vector] : prod(v) x n matrix of regimes (each row
%   represents a specific combination of states in each dimension, hence a
%   regime).
%
% - **unique_regimes** [matrix|vector] : ? x n matrix of regimes without
%   repetition and such that for each row the indexes are increasing.
%
% More About
% ------------
%
% Examples
% ---------
%
% - Regimes=mygrid(10); : one chain with 10 regimes
% - Regimes=mygrid([10,5]) : 2 chains with 5 regimes for the second one
% - Regimes=mygrid([10,5,3]): 3 chains with 3 regimes on the third one
%
% See also:

n=numel(v);
% first construct the indexes of the intervals
Regimes=[]; % first to v(1)th interval

for p=1:n
    Regimes=utils.gridfuncs.build_grid(Regimes,v(p));
end

if nargout>1
    % find the non-decreasing indices
    %--------------------------------
    % flip left right
    test2=Regimes(:,end:-1:1);
    % stamp the decreasing rows
    drow=test2(:,1:end-1)-test2(:,2:end);
    keep=~any(drow<0,2);
    unique_regimes=Regimes(keep,:);
end
