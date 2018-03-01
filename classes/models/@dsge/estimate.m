function obj=estimate(obj,varargin)
% ESTIMATE - estimates the parameters of a DSGE or RISE model
%
% More About
% ------------
%
% - ESTIMATE is the same as GENERIC/ESTIMATE except for options that
% only apply to dsge or rise models. Those options are:
%   - **estim_priors** [{[]}|struct]: This provides an alternative to
%   setting priors inside the rise/dsge model file. Each field of the
%   structure must be the name of an estimated parameter. Each field will
%   hold a cell array whose structure is described in help
%   GENERIC/setup_priors.
%
% See also: GENERIC/ESTIMATE, GENERIC/SETUP_PRIORS

if isempty(obj)
    
    mydefaults=estimate@generic(obj,varargin{:});
    
    mydefaults=[mydefaults
        {'estim_priors',[],@(x)isstruct(x),...
        'estim_priors must be a structure'}];
        
    if nargout
        
        obj=mydefaults;
        
    else
        
        clear obj
        
        disp_defaults(mydefaults);
        
    end

    
    return
    
else
    % Initially set the filtering/smoothing flag to false (during estimation).
    % This is especially important given that the objective function could be
    % optimal_simple_rule_posterior, in which case there is no filtering going
    % on.
    obj=set(obj,'kf_filtering_level',0);
    
    obj=estimate@generic(obj,varargin{:});
    
end


end
