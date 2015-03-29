function ppdata=plot_priors_and_posteriors(obj,simulation_folder)
% plot_priors_and_posteriors -- computes posterior and prior densities for
% estimated parameters 
%
% Syntax
% -------
% ::
%
%   ppdata=plot_priors_and_posteriors(obj)
%
%   ppdata=plot_priors_and_posteriors(obj,simulation_folder)
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
% Outputs
% --------
%
% - **ppdata** [struct]: optional output argument, ppdata is a structure
% containing the information needed to plot the posterior and prior
% densities. The user can always plot those using
% utils.plot.prior_posterior(ppdata.(pname)), where pname is the name of
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
    ppdata=struct();
    return
end

if nargin<2
    simulation_folder=obj.folders_paths.simulations;
end

% do posterior densitities
%---------------------------
post_dens=plot_posteriors(obj,simulation_folder);

% do prior densities for all parameters
%----------------------------------------
prior_dens=plot_priors(obj);
vnames=fieldnames(prior_dens);

% create the data
%----------------
npar=numel(vnames);
ppdata_=struct();
for ipar=1:npar
    ppdata_.(vnames{ipar})=do_one_post_prior(prior_dens.(vnames{ipar}),post_dens.(vnames{ipar}));
end

if nargout
    ppdata=ppdata_;
else
    % plot the data
    %--------------
    r0=obj.options.graphics(1);
    c0=obj.options.graphics(2);
    titel='priors and posterior marginal densities';
    
    utils.plot.multiple(@(xname)plotfunc(xname,ppdata_),...
        vnames,titel,r0,c0,...
        'FontSize',11,'FontWeight','normal');
end

    function ss=do_one_post_prior(prior,post)
        post=rmfield(post,{'x_min','x_max','tex_name'});
        
        ss=utils.miscellaneous.mergestructures(prior,post);

        % give the prior density the same range as ss.f_kdens
        %-----------------------------------------------------
        if max(ss.f_prior)==min(ss.f_prior)
            ratio=.5;
        else
            ratio=(ss.f_prior-min(ss.f_prior))/(max(ss.f_prior)-min(ss.f_prior));
        end
        ss.f_prior=min(ss.f_kdens)+ratio*(max(ss.f_kdens)-min(ss.f_kdens));
    end

end

function [tex_name,legend_]=plotfunc(pname,ppdata)
% the caller may use the tex_name information to override the title...
[~,legend_,tex_name]=utils.plot.prior_posterior(ppdata.(pname),'LineWidth',2.5);
end

