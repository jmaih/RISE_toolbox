function [pdata,hdl]=plot_posteriors(obj,simfold,parlist,...
    npoints,subset,varargin)
% Computes posterior densities for estimated parameters
%
% ::
%
%   pdata=plot_posteriors(obj)
%
%   pdata=plot_posteriors(obj,simfold)
%
%   pdata=plot_posteriors(obj,simfold,parlist)
%
%   pdata=plot_posteriors(obj,simfold,parlist,npoints)
%
%   pdata=plot_posteriors(obj,simfold,parlist,npoints,subset)
%
%   pdata=plot_posteriors(obj,simfold,parlist,npoints,subset,varargin)
%
% Args:
%
%    obj (rise | dsge | rfvar | svar): model object
%
%    simfold (empty | char | struct): location of the simulations. If
%      empty, it is assumed that the simulations are saved to disc and are
%      located in the address found in obj.folders_paths.simulations. If it is a
%      "char", this corresponds to the location of the simulation. Otherwise, if
%      it is a struct, then it has to be the output of posterior_simulator.m
%
%    npoints (numeric | {20^2}): the number of points in the
%      discretization of the prior support
%
%    parlist (empty | char | cellstr): list of the parameters for which one
%      wants to plot the posteriors
%
%    subset (cell array|{empty}): When not empty, subset is a
%     1 x 2 cell array in which the first cell contains a vector selecting
%     the columns to retain in each chain and the second column contains
%     the chains retained. Any or both of those cell array containts can be
%     empty. Whenever an entry is empty, all the information available is
%     selected. E.g. subsetting with dropping and trimming
%     mysubs={a:b:c,[1,3,5]}. In this example, the first
%     element selected is the one in position "a" and
%     thereafter every "b" element is selected until we reach
%     element in position "c". At the same time, we select
%     markov chains 1,3 and 5.
%
% Returns:
%    :
%
%    - **pdata** [struct]: optional output argument, pdata is a structure
%      containing the information needed to plot the posterior and prior
%      densities. The user can always plot those using
%      utils.plot.prior_posterior(pdata.(pname)), where pname is the name of
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

    pdata=cell(0,4);

    return

end

if nargin<5

    subset=[];

if nargin<4

    npoints=[];

    if nargin<3

        parlist=[];

        if nargin<2

            simfold=[];

        end

    end

end

end

nout=nargout;

if nout

    pdata=0;

    hdl=[];

end

nobj=numel(obj);

if nobj>1

    tmpdata=cell(1,nobj);

    tmphdl=cell(1,nobj);

    for iobj=1:nobj

        if nout

            [argouts{1:nout}]=plot_posteriors(obj(iobj),simfold,parlist,...
                npoints,subset,varargin{:});

            tmpdata{iobj}=argouts{1};

            if nout>1

                tmphdl{iobj}=argouts{2};

            end

        else

            plot_posteriors(obj(iobj),simfold,parlist,npoints,subset,varargin{:});

        end

    end

    if nout

        pdata=tmpdata;

        if nout>1

            hdl=tmphdl;

        end

    end

    return

end

% subset=[];
% 
% if iscell(simfold)
%     
%     subset=simfold{2};
%     
%     simfold=simfold{1};
%     
% end

if isempty(simfold)

    simfold=obj.folders_paths.simulations;

end

[d,~,~,~,~,~,~,best]=mcmc.reload_draws(simfold,subset);

% suppress the chain dimension

d=d(:,:);

if isempty(npoints),npoints=20^2; end

numscal=@(x)isnumeric(x) && isscalar(x);

num_fin=@(x)numscal(x) && isfinite(x);

num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;

if ~num_fin_int(npoints)

    error('npoints must be a finite and positive integer')

end

% do prior densities for all parameters
%----------------------------------------
prior_dens=plot_priors(obj,parlist);

vnames=fieldnames(prior_dens);

% locate those names
allpnames=cellfun(@(x)parser.param_texname_to_param_name(x),...
    {obj.estimation.priors.name},'uniformOutput',false);

vlocs=locate_variables(vnames,allpnames);

N=numel(prior_dens.(vnames{1}).x_prior);

is_posterior_max=isfield(obj.estimation.posterior_maximization,'mode') && ...
    ~isempty(obj.estimation.posterior_maximization.mode);

if is_posterior_max

    post_mode=obj.estimation.posterior_maximization.mode(vlocs);

end

post_mode_sim=best{1}.x(vlocs);

ff=best{1}.f;

for ib=2:numel(best)
    
    if best{ib}.f<ff
        
        ff=best{ib}.f;
        
        post_mode_sim=best{ib}.x(vlocs);
    
    end
    
end

% create the data
%----------------
npar=numel(vnames);

pdata_=struct();

% potential candidate for parallelization
for ipar=1:npar

    all_vals=d(vlocs(ipar),:);

    tex_name=prior_dens.(vnames{ipar}).tex_name;

    pdata_.(vnames{ipar})=do_one_post(ipar);

end

vargs=utils.plot.expand_varargin([],varargin{:});

if nout==0||nout==2
    % plot the data
    %--------------
    r0=obj.options.graphics(1);

    c0=obj.options.graphics(2);

    titel='posterior marginal densities';

    tmphdl=utils.plot.multiple(@(xname)plotfunc(xname,pdata_),...
        vnames,titel,r0,c0,...
        'FontSize',11,'FontWeight','normal');

end

if nout

    pdata=pdata_;

    if nout>1

        hdl=tmphdl;

    end

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

    function [tex_name,legend_]=plotfunc(pname,pdata)
        % the caller may use the tex_name information to override the title...
        [~,legend_,tex_name]=utils.plot.prior_posterior(pdata.(pname),vargs{:});

    end

end

