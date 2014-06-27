function [lnprior,retcode]=log_prior_density(obj,param)
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
        for kk=1:numel(obj.estim_distributions)
            loc=obj.estim_distrib_locations{kk};
            ld=obj.estim_distributions{kk}(param(loc),obj.estim_hyperparams(loc,1),obj.estim_hyperparams(loc,2));
            if obj.options.debug
                disp(functions(obj.estim_distributions{kk}))
                disp([{'param val','a','b','log-prior'}
                    num2cell([param(loc),obj.estim_hyperparams(loc,1),obj.estim_hyperparams(loc,2),ld])])
            end
            lnprior=lnprior+sum(ld);
            if (isnan(lnprior)||isinf(lnprior)||~isreal(lnprior))
                retcode=307;
                break
            end
        end
    elseif ischar(param) % This should be moved elsewhere to avoid recomputing the covariance...
        if ~isempty(obj.estimation.endogenous_priors) && obj.markov_chains.regimes_number==1 % not provided for switching models yet
            % this assumes the new parameter vector is already pushed in
            [tv,retcode]=theoretical_autocovariances(obj,0);
            if ~retcode
                this=obj.estimation.endogenous_priors;
                lnprior=log_density(this,tv(:,:,1));
            end
        end
    end
end

end
