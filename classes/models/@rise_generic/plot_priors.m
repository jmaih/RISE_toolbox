function pdata=plot_priors(obj,parlist,varargin)
% plot_priors -- computes prior densities for estimated parameters
%
% Syntax
% -------
% ::
%
%   ppdata=plot_priors(obj)
%
%   ppdata=plot_priors(obj,parlist)
%
%   ppdata=plot_priors(obj,parlist,varargin)
%
% Inputs
% -------
%
% - **obj** [rise|dsge|rfvar|svar]: model object
%
% - **parlist** [[]|char|cellstr]: list of the parameters for which one
% wants to plot the priors
%
% - **varargin** [pairwise plotting arguments]: 
%
% Outputs
% --------
%
% - **pdata** [struct]: optional output argument, pdata is a structure
% containing the information needed to plot the prior densities. The user
% can always plot those using utils.plot.prior_posterior(ppdata.(pname)),
% where pname is the name of one particular parameter of interest.
%
% More About
% ------------
%
% - if there are no output arguments, figures with prior densities are
% plotted, but not saved!!!.
%
% Examples
% ---------
%
% See also: utils.plot.prior_posterior


if isempty(obj)
    % For the computation of check plots, priors and posteriors
    pdata=struct('prior_discretize',20);
    return
end
if nargout
    pdata=0;
end

nobj=numel(obj);

if nobj>1
    retcode_=cell(1,nobj);
    for iobj=1:nobj
        if nargout
            retcode_{iobj}=plot_priors(obj(iobj));
        else
            plot_priors(obj(iobj));
        end
    end
    if nargout
        pdata=recode_;
    end
    return
end

%----------------------------------------
N=obj.options.prior_discretize^2;
allpnames=cellfun(@(x)parser.param_texname_to_param_name(x),...
    {obj.estimation.priors.name},'uniformOutput',false);
if nargin<2
    parlist=[];
end
if isempty(parlist)
    parlist=allpnames;
else
    if ischar(parlist)
        parlist=cellstr(parlist);
    end
    parlist=cellfun(@(x)parser.param_texname_to_param_name(x),...
        parlist,'uniformOutput',false);
end
plocs=locate_variables(parlist,allpnames);
if numel(unique(plocs))~=numel(plocs)
    error('parameter names duplicated')
end

vargs={'LineWidth',2.5};
if isempty(varargin)
    varargin=vargs;
else
    found=false;
    for ii=1:2:length(varargin)-1
        found=strcmpi(varargin{ii},vargs{1});
        if found
            break
        end
    end
    if ~found
        varargin=[vargs,varargin];
    end
end

pnames=parlist;

tex_names={obj.estimation.priors(plocs).tex_name};

distr={obj.estimation.priors(plocs).prior_distrib};
% replace the dirichlet with the beta
%-------------------------------------
distr=strrep(distr,'dirichlet','beta');
% recollect the densities
for idistr=1:numel(distr)
    distr{idistr}=distributions.(distr{idistr});
end
lb=vertcat(obj.estimation.priors(plocs).lower_bound);
ub=vertcat(obj.estimation.priors(plocs).upper_bound);
hypers=obj.estim_hyperparams(plocs,:);
npar=numel(lb);
prior_dens=struct();
for ipar=1:npar
    prior_dens.(pnames{ipar})=do_one_prior(ipar);
end
%----------------------------------------

if nargout
    pdata=prior_dens;
else
    % plot the data
    %--------------
    r0=obj.options.graphics(1);
    c0=obj.options.graphics(2);
    titel='prior densities';
    
    utils.plot.multiple(@(xname)plotfunc(xname,prior_dens),...
        pnames,titel,r0,c0,...
        'FontSize',11,'FontWeight','normal');
end

    function pdens=do_one_prior(ipar)
        pdens=struct();
        x_prior=vec(linspace(lb(ipar),ub(ipar),N));
        pdens.x_prior=x_prior;
        pdens.f_prior=distr{ipar}(x_prior,hypers(ipar,1),hypers(ipar,2));
        pdens.f_prior=exp(pdens.f_prior);
        pdens.tex_name=tex_names{ipar};
        pdens.x_min=lb(ipar);
        pdens.x_max=ub(ipar);
    end

function [tex_name,legend_]=plotfunc(pname,ppdata)
% the caller may use the tex_name information to override the title...
[~,legend_,tex_name]=utils.plot.prior_posterior(ppdata.(pname),varargin{:});
end

end
