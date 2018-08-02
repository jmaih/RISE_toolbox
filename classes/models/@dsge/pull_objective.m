function varargout=pull_objective(obj,varargin)
% INTERNAL FUNCTION: pulls the objective function to optimize for a DSGE or RISE model
%
% Note:
%
%    - pull_objective is the same as generic/pull_objective except for
%      the kf_filtering_level
%
% See also:
%    - rise_generic/pull_objective
%

if ~isempty(obj)

    obj=set(obj,'kf_filtering_level',0);

end

[varargout{1:nargout}]=pull_objective@generic(obj,varargin{:});

end