function [log_post,log_lik,log_prior,Incr,retcode,obj]=log_posterior_kernel(obj,param)
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
    log_post=struct();
    return
end
% estim_hyperparams=obj.estim_hyperparams;
log_prior=nan(1,2);
log_lik=-obj.options.estim_penalty;
log_post=log_lik;
Incr=[];

likelihood_func=obj.routines.likelihood;
[log_prior(1),retcode]=log_prior_density(obj,param);
if ~retcode
    [log_lik,Incr,retcode,obj]=likelihood_func(param,obj);
    if ~retcode
        [log_prior(2),retcode]=log_prior_density(obj,'endogenous');
        if ~retcode
            log_post=log_lik+sum(log_prior);
            if log_post<=-obj.options.estim_penalty
                retcode=306; % unlikely parameter vector
            end
        end
    end
    % [LogLik,Incr,retcode,obj]=likelihood_dsge_var(params,obj,kf_filtering_level)
    % [LogLik,Incr,retcode,obj]=likelihood_markov_switching_dsge(params,obj,kf_filtering_level)
    % [LogLik,Incr,retcode,obj]=likelihood_optimal_simple_rule(params,obj,kf_filtering_level)
end
if obj.options.debug
    disp(['log_lik ',num2str(log_lik)])
    disp(['log_prior ',num2str(log_prior(1))])
    disp(['log_endo_prior ',num2str(log_prior(2))])
    disp(['log_post ',num2str(log_post)])
    utils.error.decipher(retcode)
    if obj.estimation_under_way
        error([mfilename,':: no debugging of this function should occur during estimation, at least I do not see why'])
    else
        keyboard
    end
end

end

