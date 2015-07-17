function [dens,N,pnames]=compute_prior_densities(obj)
% compute_prior_densities -- computes prior densities over a discretized
% number of values over the range of each estimated parameter
%
% Syntax
% -------
% ::
%
%   [dens,N,pnames]=compute_prior_densities(obj)
%
% Inputs
% -------
%
% - **obj** [rise|dsge|rfvar|svar]: model object
%
% Outputs
% --------
%
% - **dens** [struct]: structure whose fields are the names of the
% estimated parameters. Those names are follows are followed by a structure
% with fields:
%   - **x_prior** [vector]: x-axis values discretized from the
%       range/support of the corresponding parameter
%   - **f_prior** [vector]: y-axis values of densities of x_prior
%
% - **N** [scalar]: number of discretized values
%
% - **pnames** [cellstr]: names of the parameters
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

N=obj.options.prior_discretize^2;

pnames=cellfun(@(x)parser.param_name_to_valid_param_name(x),...
    {obj.estimation.priors.name},'uniformOutput',false);

distr={obj.estimation.priors.prior_distrib};
% recollect the densities
for idistr=1:numel(distr)
    distr{idistr}=distributions.(distr{idistr});
end
lb=vertcat(obj.estimation.priors.lower_bound);
ub=vertcat(obj.estimation.priors.upper_bound);
hypers=obj.estim_hyperparams;
npar=numel(lb);

dens=struct();
for ipar=1:npar
    x_prior=vec(linspace(lb(ipar),ub(ipar),N));
    dens.(pnames{ipar}).x_prior=x_prior;
    dens.(pnames{ipar}).f_prior=distr{ipar}(x_prior,hypers(ipar,1),hypers(ipar,2));
end

end
