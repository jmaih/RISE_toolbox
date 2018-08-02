function [lnprior,retcode]=log_prior_density(obj, param)
% Computes the probability density function of the prior corresponding to the parameter values
%
% Note:
%    In effort to make RISE modular, this function is available so that one can use a different sampler if needed, but most likely, one should just use available stock samplers.
%
% ::
%
%    [lnprior, retcode] = log_prior_density(model, param)
%
% Args:
%    model (dsge | rise object): model object
%    param (column vector): parameter values
%
% Returns:
%    : [lnprior, retcode]
%
%    - **lnprior** (double): log of prior density function
%    - **retcode**: return code
%
% See also:
%    - :func:`log_posterior_kernel <dsge.log_posterior_kernel>`
%

nobj=numel(obj);

if nobj==0

    mydefaults=the_defaults();

    if nargout

        lnprior=mydefaults;

    else

        disp_defaults(mydefaults);

    end

    return

end

if nargin<2

    param=[];

end

lnprior=0;

retcode=0;

epdata=obj.estim_priors_data;

if ~isempty(param)

    if isnumeric(param)

        [lnprior,retcode]=utils.estim.prior_evaluation_engine(epdata,param,...
            lnprior);

    elseif ischar(param)

        if ~isempty(obj.estimation.endogenous_priors)
            % do the uncorrelated guys
            %--------------------------
            pp=obj.options.estim_endogenous_priors(obj);

            [lnprior,retcode]=utils.estim.prior.evaluate_uncorrelated(...
                obj.estim_endogenous_priors_data,pp,lnprior);

            if retcode

                retcode=309;

            end

        end

    end

end

end

function d=the_defaults()

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x);

d={
    'prior_trunc',1e-10,@(x)num_fin(x) && x>0 && x<1,'prior_trunc must be in (0,1)'
    };

end
