function [init,sampler,total_draws]=initialize_posterior_simulation(obj)
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

lb=[obj.estimation.priors.lower_bound]';
ub=[obj.estimation.priors.upper_bound]';
drawfun=@(x,cCs)truncated_multivariate_normal.quick_and_dirty(x,cCs,lb,ub);
sampler=@utils.mcmc.mh_sampler;
c=obj.options.mcmc_initial_covariance_tune;
if c<=0
    error([mfilename,':: mcmc_initial_covariance_tune must be positive'])
end
vcov=utils.cov.nearest(obj.estimation.posterior_maximization.vcov);
x0=obj.estimation.posterior_maximization.mode;
not_optimized=isempty(x0);
if not_optimized
    x0=[obj.estimation.priors.start]';
    npar=numel(obj(1).estimation.priors);
    vcov=eye(npar);
    [minus_log_post,retcode]=pull_objective(x0);
    f0=minus_log_post;
    if retcode
        msg=utils.error.decipher(retcode);
        error(['At MCMC simulation, starting values leads to the following problem ''',...
            msg,...
            ''' .This can be corrected. please contact junior.maih@gmail.com'])
    end
else
    f0=-obj.estimation.posterior_maximization.log_post;
end
fprintf(1,'%s %8.4f\n','Initial value of log posterior ',f0);
% adapt the covariance matrix automatically if we have a
% diagonal original covariance matrix regardless of ...
%--------------------------------------------------------------
adapt_covariance=all(all(abs(diag(diag(vcov))-vcov)<1e-12));
adapt_covariance=adapt_covariance||obj.options.mcmc_adapt_covariance;
%----------------------------------------
CS=transpose(chol(vcov));
%----------------------------------------
number_of_burns=round(obj.options.mcmc_burn_rate*obj.options.mcmc_number_of_simulations);
init=struct('x0',x0,'f0',f0,'CS',CS,...
    'c',c,'objective',pull_objective(obj),...
    'drawfun',@(varargin)drawfun(varargin{:}),...
    'adapt_covariance',adapt_covariance,...
    'thin',obj.options.mcmc_thin,...
    'burnin',number_of_burns,...
    'funevals',0);
total_draws=number_of_burns+obj.options.mcmc_number_of_simulations;

end
