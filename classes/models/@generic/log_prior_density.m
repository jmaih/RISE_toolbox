function [lnprior,retcode]=log_prior_density(obj,param)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

nobj=numel(obj);
if nobj==0
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    lnprior=struct('prior_trunc',1e-10,...
        'prior_endogenous',false);
    % alternatives are 
    %           []: same as false --> no endogenous priors
    %           true: endogenous priors with all observables and
    %                 strength=sample size
    %           struct('targets',target_names,'prior_sample',integer): the
    %                 listed target variables must be observable and the 
    %                 prior sample represents the strength of the beliefs
    return
end
if nargin<2
    param=[];
end
        
lnprior=0;  
%
retcode=0;
if ~isempty(param)
    if isnumeric(param)
        % do the uncorrelated guys
        %--------------------------
        for kk=1:numel(obj.estim_distributions)
            if ischar(obj.estim_distributions{kk})
                % that is the dirichlet,skip it
                continue
            end
            loc=obj.estim_distrib_locations{kk};
            a=obj.estim_hyperparams(loc,1);
            b=obj.estim_hyperparams(loc,2);
            ld=obj.estim_distributions{kk}(param(loc),a,b);
            if obj.options.debug
                disp(functions(obj.estim_distributions{kk}))
                disp([{'param val','a','b','log-prior'}
                    num2cell([param(loc),a,b,ld])])
            end
            lnprior=lnprior+sum(ld);
            if (isnan(lnprior)||isinf(lnprior)||~isreal(lnprior))
                retcode=307;
                break
            end
        end
        % do the dirichlet guys
        %-----------------------
        if ~retcode && ~isempty(obj.estim_dirichlet)
            ndirich=numel(obj.estim_dirichlet);
            for id=1:ndirich
                loc=obj.estim_dirichlet(id).location;
                lnprior=lnprior+obj.estim_dirichlet(id).lpdfn(param(loc));
                if (isnan(lnprior)||isinf(lnprior)||~isreal(lnprior))
                    retcode=307;
                    break
                end
            end
        end
    elseif ischar(param) % This should be moved elsewhere to avoid recomputing the covariance...
        if ~isempty(obj.estimation.endogenous_priors) && obj.markov_chains.regimes_number==1 % not provided for switching models yet
            [endo_priors_func,vargs]=utils.code.user_function_to_rise_function(obj.estimation.endogenous_priors);
            [lnprior,obj.user_endo_priors_info,retcode]=endo_priors_func(obj,obj.user_endo_priors_info,vargs{:});
            if retcode
                retcode=309;
            end
        end
    end
end

end
