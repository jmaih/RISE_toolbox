function [ppdata,hdl]=plot_priors_and_posteriors(obj,sim_fold,...
    parlist,trunc,npoints,varargin)
% Compute posterior and prior densities for estimated parameters
%
% ::
%
%   ppdata=plot_priors_and_posteriors(obj)
%   ppdata=plot_priors_and_posteriors(obj,sim_fold)
%   ppdata=plot_priors_and_posteriors(obj,sim_fold,parlist)
%   ppdata=plot_priors_and_posteriors(obj,sim_fold,parlist,trunc)
%   ppdata=plot_priors_and_posteriors(obj,sim_fold,parlist,trunc,npoints)
%   ppdata=plot_priors_and_posteriors(obj,sim_fold,parlist,trunc,npoints,varargin)
%
% Args:
%
%    obj (rise | dsge | rfvar | svar): model object
%
%    sim_fold (empty | char | struct): location of the simulations. If
%      empty, it is assumed that the simulations are saved to disc and are
%      located in the address found in obj.folders_paths.simulations. If it is a
%      "char", this corresponds to the location of the simulation. Otherwise, if
%      it is a struct, then it has to be the output of posterior_simulator.m
%
%    parlist (empty | char | cellstr): list of the parameters for which one
%      wants to plot the priors and the posteriors
%
%    trunc (numeric | {1e-3}): serves to truncate the support
%
%    npoints (numeric | {20^2}): the number of points in the
%      discretization of the prior support
%
% Returns:
%    :
%
%    - **ppdata** [struct]: optional output argument, ppdata is a structure
%      containing the information needed to plot the posterior and prior
%      densities. The user can always plot those using
%      utils.plot.prior_posterior(ppdata.(pname)), where pname is the name of
%      one particular parameter of interest.
%
%    - **hdl** [vector]: optional output argument, placeholder for handles to
%      graphs
%
% Note:
%
%    - if there are no output arguments, figures with posterior and prior
%      marginal densities are plotted, but not saved!!!.
%      see also utils.plot.prior_posterior
%

if isempty(obj)

    ppdata=cell(0,4);

    return

end

    if nargin<5

        npoints=[];

        if nargin<4

            trunc=[];

            if nargin<3

                parlist=[];

                if nargin<2

                    sim_fold=[];

                end

            end

        end

    end

nout=nargout;

if nout

    ppdata=0;

    hdl=[];

end

nobj=numel(obj);

if nobj>1

    tmpdata=cell(1,nobj);

    tmphdl=cell(1,nobj);

    for iobj=1:nobj

        if nout

            [argouts{1:nout}]=plot_priors_and_posteriors(obj(iobj),...
                sim_fold,parlist,trunc,npoints,varargin{:});

            tmpdata{iobj}=argouts{1};

            if nout>1

                tmphdl{iobj}=argouts{2};

            end

        else

            plot_priors_and_posteriors(obj(iobj),sim_fold,...
                parlist,trunc,npoints,varargin{:});

        end

    end

    if nout

        ppdata=tmpdata;

        if nout>1

            hdl=tmphdl;

        end

    end

    return

end

if isempty(sim_fold)

    sim_fold=obj.folders_paths.simulations;

end

% do posterior densities
%---------------------------
post_dens=plot_posteriors(obj,sim_fold,parlist,npoints);

vnames=fieldnames(post_dens);

% do prior densities for all parameters
%----------------------------------------
prior_dens=plot_priors(obj,vnames,trunc,npoints);

% create the data
%----------------
npar=numel(vnames);

ppdata_=struct();

for ipar=1:npar

    ppdata_.(vnames{ipar})=do_one_post_prior(prior_dens.(vnames{ipar}),...
        post_dens.(vnames{ipar}));

end

vargs=utils.plot.expand_varargin([],varargin{:});

if nout==0||nout==2
    % plot the data
    %--------------
    r0=obj.options.graphics(1);

    c0=obj.options.graphics(2);

    titel='prior and posterior(marginal) densities';

    tmphdl=utils.plot.multiple(@(xname)plotfunc(xname,ppdata_),...
        vnames,titel,r0,c0,...
        'FontSize',11,'FontWeight','normal');

end

if nout

    ppdata=ppdata_;

    if nout>1

        hdl=tmphdl;

    end

end

    function ss=do_one_post_prior(prior,post)

        prior.x_min=min(prior.x_min,post.x_min);

        prior.x_max=max(prior.x_max,post.x_max);

        post=rmfield(post,{'x_min','x_max','tex_name'});

        ss=utils.miscellaneous.mergestructures(prior,post);

    end

    function [tex_name,legend_]=plotfunc(pname,ppdata)
        % the caller may use the tex_name information to override the title...
        [~,legend_,tex_name]=utils.plot.prior_posterior(ppdata.(pname),vargs{:});
    end

end

