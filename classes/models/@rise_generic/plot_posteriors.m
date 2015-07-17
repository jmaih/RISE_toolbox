function pdata=plot_posteriors(obj,simulation_folder,parlist)
% plot_posteriors -- computes posterior densities for estimated parameters
%
% Syntax
% -------
% ::
%
%   pdata=plot_posteriors(obj)
%
%   pdata=plot_posteriors(obj,simulation_folder)
%
% Inputs
% -------
%
% - **obj** [rise|dsge|rfvar|svar]: model object
%
% - **simulation_folder** [[]|char|struct]: location of the simulations. If
% empty, it is assumed that the simulations are saved to disc and are
% located in the address found in obj.folders_paths.simulations. If it is a
% "char", this corresponds to the location of the simulation. Otherwise, if
% it is a struct, then it has to be the output of posterior_simulator.m
%
% - **parlist** [[]|char|cellstr]: list of the parameters for which one
% wants to plot the posteriors
%
% Outputs
% --------
%
% - **pdata** [struct]: optional output argument, pdata is a structure
% containing the information needed to plot the posterior and prior
% densities. The user can always plot those using
% utils.plot.prior_posterior(pdata.(pname)), where pname is the name of
% one particular parameter of interest. 
%
% More About
% ------------
%
% - if there are no output arguments, figures with posterior and prior
% marginal densities are plotted, but not saved!!!.
% see also utils.plot.prior_posterior
%
% Examples
% ---------
%
% See also: 

if isempty(obj)
    pdata=struct();
    return
end

if nargin<3
    parlist=[];
    if nargin<2
        simulation_folder=[];
    end
end

if isempty(simulation_folder)
    simulation_folder=obj.folders_paths.simulations;
end

is_saved_to_disk=ischar(simulation_folder);
if is_saved_to_disk
    W = what(simulation_folder);
    W=W.mat;
    locs=find(strncmp('chain_',W,6));
    if isempty(locs)
        error([mfilename,':: no simulations found'])
    end
    W=strrep(W(locs),'.mat','');
elseif isstruct(simulation_folder)
    W=fieldnames(simulation_folder);
else
    error('wrong specification of input')
end
number_of_matrices=numel(W);

% do prior densities for all parameters
%----------------------------------------
prior_dens=plot_priors(obj,parlist);
vnames=fieldnames(prior_dens);
% locate those names
allpnames=cellfun(@(x)parser.param_name_to_valid_param_name(x),...
    {obj.estimation.priors.name},'uniformOutput',false);
vlocs=locate_variables(vnames,allpnames);
N=numel(prior_dens.(vnames{1}).x_prior);

is_posterior_max=isfield(obj.estimation.posterior_maximization,'mode');
post_mode_sim=[];
f_post_mode_sim=-inf;
if is_posterior_max
    post_mode=obj.estimation.posterior_maximization.mode(vlocs);
%     f_post_mode=obj.estimation.posterior_maximization.log_post;
end

% create the data
%----------------
npar=numel(vnames);
pdata_=struct();
for ipar=1:npar
    all_vals=[];
    for m=1:number_of_matrices
        if is_saved_to_disk
            tmp=load([simulation_folder,filesep,W{m}]);
        else
            tmp=simulation_folder.(W{m});
        end
        if ipar==1
            % try and locate the sampling posterior mode
            fm=-tmp.minus_logpost_params;
            best=find(fm==max(fm),1,'first');
            if fm(best)>f_post_mode_sim
                post_mode_sim=tmp.Params(vlocs,best);
                f_post_mode_sim=fm(best);
            end
        end
        all_vals=[all_vals;tmp.Params(vlocs(ipar),:).']; %#ok<AGROW>
    end
    tex_name=prior_dens.(vnames{ipar}).tex_name;
    pdata_.(vnames{ipar})=do_one_post(ipar);
end

if nargout
    pdata=pdata_;
else
    % plot the data
    %--------------
    r0=obj.options.graphics(1);
    c0=obj.options.graphics(2);
    titel='posterior marginal densities';
    
    utils.plot.multiple(@(xname)plotfunc(xname,pdata_),...
        vnames,titel,r0,c0,...
        'FontSize',11,'FontWeight','normal');
end

    function ss=do_one_post(ipar)
        ss=struct();
        ss.mean_sim=mean(all_vals);
        x_min_sim = min(all_vals);
        x_max_sim = max(all_vals);
        [ss.f_kdens,ss.x_kdens]=distributions.kernel_density(all_vals,[],[],'normal',N);
        if is_posterior_max
            ss.post_mode=post_mode(ipar);
        end
        ss.post_mode_sim=post_mode_sim(ipar);
        ss.x_min=x_min_sim;
        ss.x_max=x_max_sim;
        ss.tex_name=tex_name;
    end

end

function [tex_name,legend_]=plotfunc(pname,pdata)
% the caller may use the tex_name information to override the title...
[~,legend_,tex_name]=utils.plot.prior_posterior(pdata.(pname),'LineWidth',2.5);
end

