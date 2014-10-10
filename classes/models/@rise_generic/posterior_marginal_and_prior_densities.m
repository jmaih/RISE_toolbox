function ppdata=posterior_marginal_and_prior_densities(obj,simulation_folder)
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
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

% computes the posterior and marginal densities.
% - the optional argument simulation_folder can either be the output of
% posterior_simulator.m or the path where the posterior simulations are
% stored.
% - the optional output argument ppdata is a structure containing the
% information needed to plot the posterior and prior densities. The user
% can always plot those using utils.plot.prior_posterior(ppdata.(pname))
% where pname is the name of one particular parameter of interest.
% - if there are no output arguments, figures with posterior and prior
% marginal densities are plotted, but not saved!!!.
% see also utils.plot.prior_posterior
if isempty(obj)
    ppdata=struct();
    return
end

if nargin<2
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


N=obj.options.prior_discretize^2;

distr={obj.estimation.priors.prior_distrib};
% recollect the densities
for idistr=1:numel(distr)
    distr{idistr}=distributions.(distr{idistr});
end
lb=vertcat(obj.estimation.priors.lower_bound);
ub=vertcat(obj.estimation.priors.upper_bound);

post_mode=obj.estimation.posterior_maximization.mode;
f_post_mode=obj.estimation.posterior_maximization.log_post;
post_mode_sim=post_mode;
f_post_mode_sim=f_post_mode;

hypers=obj.estim_hyperparams;

vnames=cellfun(@(x)parser.param_name_to_valid_param_name(x),...
    {obj.estimation.priors.name},'uniformOutput',false);
tex_names={obj.estimation.priors.tex_name};

% create the data
%----------------
npar=size(post_mode,1);
ppdata_=struct();
for ipar=1:npar
    all_vals=[];
    for m=1:number_of_matrices
        if is_saved_to_disk
            tmp=load([simulation_folder,filesep,W{m}]);
        else
            tmp=simulation_folder.(W{m});
        end
        Params=tmp.Params(ipar,:);
        if ipar==1
            % try and locate the sampling posterior mode
            fm=-tmp.minus_logpost_params;
            best=find(fm==max(fm),1,'first');
            if fm(best)>f_post_mode_sim
                post_mode_sim=Params(:,best);
                f_post_mode_sim=fm(best);
            end
        end
        all_vals=[all_vals;Params(:)]; %#ok<AGROW>
    end
    ppdata_.(vnames{ipar})=do_one_post_prior(ipar);
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

    function ss=do_one_post_prior(ipar)
        ss=struct();
        ss.tex_name=tex_names{ipar};
        ss.mean_sim=mean(all_vals);
        ss.min_sim = min(all_vals);
        ss.max_sim = max(all_vals);
        [ss.f_kdens,ss.x_kdens]=distributions.kernel_density(all_vals,[],[],'normal',N);
        ss.post_mode=post_mode(ipar);
        ss.post_mode_sim=post_mode_sim(ipar);
        ss.x_prior=linspace(lb(ipar),ub(ipar),N);
        ss.x_prior=ss.x_prior(:);
        ss.f_prior=distr{ipar}(ss.x_prior,hypers(ipar,1),hypers(ipar,2));
        
        % give it the same range as ss.f_kdens
        %-------------------------------------
        if max(ss.f_prior)==min(ss.f_prior)
            ratio=.5;
        else
            ratio=(ss.f_prior-min(ss.f_prior))/(max(ss.f_prior)-min(ss.f_prior));
        end
        ss.f_prior=min(ss.f_kdens)+ratio*(max(ss.f_kdens)-min(ss.f_kdens));
    end

end

function [tex_name,legend_]=plotfunc(pname,ppdata)
tex_name=ppdata.(pname).tex_name;
[~,legend_]=utils.plot.prior_posterior(ppdata.(pname),'LineWidth',2.5);
end

