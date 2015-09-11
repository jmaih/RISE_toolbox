function varargout=pull_objective(obj,varargin)
% PULL_OBJECTIVE -- pulls the objective function to optimize
%
% More About
% ------------
%
% - PULL_OBJECTIVE is the same as RISE_GENERIC/PULL_OBJECTIVE except for
% the fact that constant-parameter RFVAR models sometimes have analytical
% solutions. In that case the posterior has a known solution. Hence the
% simulation exercise is greatly simplified. In the constant-parameter
% case then, option **vp_analytical_post_mode** can be set to true or false
% with the former being the default.
%
% See also: RISE_GENERIC/PULL_OBJECTIVE, DSGE/PULL_OBJECTIVE

nout=nargout;
if isempty(obj)
    [varargout{1:nout}]=pull_objective@svar(obj,varargin{:});
    % use analytical posterior whenever possible. Normally this should have
    % been used in some estimation method but I do not want to write a
    % specialized one just for the sake of this one element.
    varargout{1}.vp_analytical_post_mode=true;
    return
end

if ~isempty(varargin)
    obj=set(obj,varargin{:});
end

if obj.markov_chains.regimes_number==1 && obj.options.vp_analytical_post_mode
    % if the prior is set, the following should available even if the model has
    % not been estimated yet.
    lb=[obj.estimation.priors.lower_bound]';
    ub=[obj.estimation.priors.upper_bound]';
    x0=obj.estimation.posterior_maximization.mode;
    vcov=obj.estimation.posterior_maximization.vcov;
    mycell={obj.constant_var_data,lb,ub,x0,vcov};
    varargout=mycell(1:nout);
else
    [varargout{1:nout}]=pull_objective@svar(obj,varargin{:});
end

end