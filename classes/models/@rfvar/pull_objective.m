function varargout=pull_objective(obj,varargin)
% PULL_OBJECTIVE -- pulls the objective function to optimize
%
% Syntax
% -------
% ::
%
%   [ff,lb,ub]=PULL_OBJECTIVE(obj)
%
%   [ff,lb,ub]=PULL_OBJECTIVE(obj,varargin)
%
%   [ff,lb,ub,x0]=PULL_OBJECTIVE(obj,varargin)
%
%   [ff,lb,ub,x0,vcov]=PULL_OBJECTIVE(obj,varargin)
%
% Inputs
% -------
%
% - **obj** [rise|dsge|svar|rfvar]: initial model object
%
% - **varargin** [pairwise addional inputs]: usual RISE arguments
%
% Outputs
% --------
%
% - **ff** [function handle]: for computing "minus log posterior kernel"
%
% - **lb** [d x 1 vector]: lower bound of the parameters to optimize
%
% - **ub** [d x 1 vector]: upper bound of the parameters to optimize
%
% - **x0** [d x 1 vector]: posterior mode if available
%
% - **vcov** [d x d matrix]: covariance matrix at the posterior mode if
% available.
%
% More About
% ------------
%
% - The function can be used for :
%   - optimization,
%   - gradient computation,
%   - hessian computation,
%   - posterior simulation
%
% - Using this function is potentially costly, one could alternatively
% simply use log_posterior_kernel. However, if there are restrictions, they
% will not be enforced. Nevertheless it is an interesting proposition that
% should be investigated further.
%
% Examples
% ---------
%
% See also:

nout=nargout;
if isempty(obj)
    varargout=cell(1,nout);
    varargout{1}=struct();
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