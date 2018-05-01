function [Regimes,unique_regimes]=mygrid(v,vectorized)
% mygrid creates a grid of points
%
% ::
%
%
%   [Regimes,unique_regimes]=mygrid(v)
%
%   [Regimes,unique_regimes]=mygrid(v,vectorized)
%
% Args:
%
%    - **v** [scalar|vector] : n x 1 or 1 x n each entry in the vector
%      represents the number of elements along that dimension
%
%    - **vectorized** [{true}|false] : if true, use a vectorization instead of
%    a for loop
%
% Returns:
%    :
%
%    - **Regimes** [matrix|vector] : prod(v) x n matrix of regimes (each row
%      represents a specific combination of states in each dimension, hence a
%      regime).
%
%    - **unique_regimes** [matrix|vector] : ? x n matrix of regimes without
%      repetition and such that for each row the indexes are increasing.
%
% Note:
%
% Example:
%
%    - Regimes=mygrid(10); : one chain with 10 regimes
%    - Regimes=mygrid([10,5]) : 2 chains with 5 regimes for the second one
%    - Regimes=mygrid([10,5,3]): 3 chains with 3 regimes on the third one
%
%    See also:
if nargin<2
    
    vectorized=true;
    
end

n=numel(v);
% first construct the indexes of the intervals
Regimes=[]; % first to v(1)th interval

for p=1:n
    
    Regimes=utils.gridfuncs.build_grid(Regimes,v(p),vectorized);
    
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
