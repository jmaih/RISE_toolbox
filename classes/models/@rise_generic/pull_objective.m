function [ff,lb,ub,x0,vcov,obj]=pull_objective(obj,varargin)
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
%   [ff,lb,ub,x0,vcov,obj]=PULL_OBJECTIVE(obj,varargin)
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
% - **obj** [rise|dsge|svar|rfvar]: updated model object
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
% - The updated object should be used for doing various exercises (irfs,
% simulations, etc.) if the posterior mode is not computed.
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

if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    ff=struct();
    return
end

if ~isempty(varargin)
    obj=set(obj,varargin{:});
end
% if the prior is set, the following should available even if the model has
% not been estimated yet.
lb=[obj.estimation.priors.lower_bound]';
ub=[obj.estimation.priors.upper_bound]';
x0=obj.estimation.posterior_maximization.mode;
vcov=obj.estimation.posterior_maximization.vcov;

% same as in estimate lines 266-267
%-----------------------------------
obj=setup_linear_restrictions(obj);
obj=setup_general_restrictions(obj);

% shorten everything in the presence of linear restrictions
%-----------------------------------------------------------
[x0,vcov,lbub_short]=shorten_under_linear_restrictions(obj,x0,{vcov},[lb,ub]);

lb=lbub_short(:,1);
ub=lbub_short(:,2);

% try to avoid unnecessary computations like storing filters and so on
%----------------------------------------------------------------------
obj.estimation_under_way=true;
obj.options.kf_filtering_level=0;

ff=@engine;

    function [minus_log_post,retcode]=engine(x)
        retcode=0;
        if any(x<lb)||any(x>ub)
            retcode=7;
        else
            [~,minus_log_post,~,issue,viol]=estimation_wrapper(obj,[],x,lb,ub,0);
            % function evaluations computed elsewhere
            if ~isempty(issue)
                retcode=issue;
            end
            if ~isempty(viol) && any(viol>0)
                retcode=7;
            end
        end
        if retcode
            minus_log_post=obj.options.estim_penalty;
        end
    end
end