function [pdata,hdl]=plot_posteriors(obj,simulation_folder,parlist,varargin)
% plot_posteriors -- computes posterior densities for estimated parameters
%
% ::
%
%
%   pdata=plot_posteriors(obj)
%
%   pdata=plot_posteriors(obj,simulation_folder)
%
% Args:
%
%    - **obj** [rise|dsge|rfvar|svar]: model object
%
%    - **simulation_folder** [empty|char|struct]: location of the simulations. If
%    empty, it is assumed that the simulations are saved to disc and are
%    located in the address found in obj.folders_paths.simulations. If it is a
%    "char", this corresponds to the location of the simulation. Otherwise, if
%    it is a struct, then it has to be the output of posterior_simulator.m
%
%    - **parlist** [empty|char|cellstr]: list of the parameters for which one
%    wants to plot the posteriors
%
% Returns:
%    :
%
%    - **pdata** [struct]: optional output argument, pdata is a structure
%    containing the information needed to plot the posterior and prior
%    densities. The user can always plot those using
%    utils.plot.prior_posterior(pdata.(pname)), where pname is the name of
%    one particular parameter of interest.
%
%    - **hdl** [vector]: optional output argument, placeholder for handles to
%    graphs
%
% Note:
%
%    - if there are no output arguments, figures with posterior and prior
%    marginal densities are plotted, but not saved!!!.
%    see also utils.plot.prior_posterior
%
% Example:
%
%    See also:

if isempty(obj)
    
    pdata=cell(0,4);
    
    return
    
end

if nargin<3
    
    parlist=[];
    
    if nargin<2
        
        simulation_folder=[];
        
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
            
            [argouts{1:nout}]=plot_posteriors(obj(iobj),simulation_folder,parlist,varargin{:});
            
            tmpdata{iobj}=argouts{1};
            
            if nout>1
                
                tmphdl{iobj}=argouts{2};
                
            end
            
        else
            
            plot_posteriors(obj(iobj),simulation_folder,parlist,varargin{:});
            
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


if isempty(simulation_folder)
    
    simulation_folder=obj.folders_paths.simulations;
    
end

is_saved_to_disk=ischar(simulation_folder);

number_of_matrices=1;

if is_saved_to_disk
    
    W = what(simulation_folder);
    
    W=W.mat;
    
    W=strrep(W(locs),'.mat','');
    
    number_of_matrices=numel(W);
    
elseif ~isstruct(simulation_folder)
    
    error('wrong specification of simulation_folder')
    
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

post_mode_sim=[];

f_post_mode_sim=-inf;

if is_posterior_max
    
    post_mode=obj.estimation.posterior_maximization.mode(vlocs); 
    
end

% create the data
%----------------
npar=numel(vnames);

pdata_=struct();

% potential candidate for parallelization
for ipar=1:npar
    
    all_vals=[];
    
    for m=1:number_of_matrices
        
        if is_saved_to_disk
            
            tmp=load([simulation_folder,filesep,W{m}]);
            
        else
            
            tmp=simulation_folder;
            
        end
        
        if ipar==1
            % try and locate the sampling posterior mode
            fm=-[tmp.f];
            
            best=find(fm==max(fm),1,'first');
            
            if fm(best)>f_post_mode_sim
                
                post_mode_sim=tmp(best).x(vlocs);
                
                f_post_mode_sim=fm(best);
                
            end
            
        end
        
		Params=[tmp.x];
        
        all_vals=[all_vals;Params(vlocs(ipar),:).']; %#ok<AGROW>
        
    end
    
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

