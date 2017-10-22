function varargout=pull_objective(obj,varargin)
% PULL_OBJECTIVE -- pulls the objective function to optimize for a DSGE or
% RISE model
%
% More About
% ------------
%
% - PULL_OBJECTIVE is the same as GENERIC_SWITCH/PULL_OBJECTIVE except for
% the kf_filtering_level
%
% Examples
% ---------
%
% See also: RISE_GENERIC/PULL_OBJECTIVE

if ~isempty(obj)
    
    obj=set(obj,'kf_filtering_level',0);
    
end

[varargout{1:nargout}]=pull_objective@generic(obj,varargin{:});

end