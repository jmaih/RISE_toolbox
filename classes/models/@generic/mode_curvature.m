function db=mode_curvature(obj,varlist,type)
% mode_curvature -- checks the curvature at the posterior mode
%
% ::
%
%
%   db=mode_curvature(obj)
%
%   db=mode_curvature(obj,varlist)
%
%   db=mode_curvature(obj,varlist,type)
%
% Args:
%
%    - **obj** [rise|dsge|rfvar|svar]: model object
%
%    - **varlist** [char|cellstr|empty]: list of parameters for which we want
%    to check the curvature
%
%    - **type** [{'max'}|'min'|'range']: normalization of the log-posterior
%    and the log-likelihood.
%
% Returns:
%    :
%
%    - **db** [struct]: structure containing the information to plot the
%    curvature. Each field is the name of a particular parameter
%
% Note:
%
%    - when no output is requested, plots are made but not saved.
%
%    - one way to plot the curvatures from the output is to use the function
%    utils.plot.curvature
%
% Example:
%
%    See also: utils.plot.curvature

if isempty(obj)
    
    if nargout>1
        
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
        
    end
    
    db=cell(0,4);
    
    return
    
end

if nargin<3
    
    type=[];
    
    if nargin<2
        
        varlist=[];
        
    end
    
end

if isempty(type)
    
    type='max';
    
end

AllNames={obj.estimation.priors.name};

if isempty(varlist)
    
    varlist=AllNames;
    
end

if ischar(varlist)
    
    varlist=cellstr(varlist);
    
end

list_locs=locate_variables(varlist,AllNames);

N=obj.options.prior_discretize;

LB=vertcat(obj.estimation.priors.lower_bound);

UB=vertcat(obj.estimation.priors.upper_bound);

xmode=obj.estimation.posterior_maximization.mode;

obj.options.kf_filtering_level=0; % do not filter

npar=numel(list_locs);

db_=cell(2,npar);

mainfunc=@do_one_parameter;

varnames={obj.estimation.priors(list_locs).name};

for ipar=1:npar
    
    var_id=list_locs(ipar);
    
    vname=varnames{ipar};
    
    db_(:,ipar)={vname,mainfunc(var_id)}';
    
end

fields=db_(1,:);

db_=cell2struct(db_(2,:),fields,2);

plotit=nargout==0;

if plotit
    % plot the data
    %--------------
    r0=obj.options.graphics(1);
    
    c0=obj.options.graphics(2);
    
    titel='Curvature around the found mode';
    
    utils.plot.multiple(@(xname)do_one_plot(xname,db_),...
        varnames,titel,r0,c0,...
        'FontSize',11,'FontWeight','normal');
    
else
    
    db=db_;
    
end

    function pp=do_one_parameter(var_id)
        
        pp=struct();
        
        vtexname=obj.estimation.priors(var_id).tex_name;
        
        pp.tex_name=vtexname;
        
        pp.mode=xmode(var_id);
        
        pp.log_post_mode=obj.estimation.posterior_maximization.log_post;
        
        pp.log_lik_mode=obj.estimation.posterior_maximization.log_lik;
        
        low = max(LB(var_id),0.8*xmode(var_id));
        
        high = min(UB(var_id),1.2*xmode(var_id));
        
        pp.x=sort([linspace(low,high,N),xmode(var_id)]);
        
        posj=find(abs(pp.x-xmode(var_id))==min(abs(pp.x-xmode(var_id))),1,'first');
        
        pp.log_post=zeros(1,N+1);
        
        pp.log_lik=zeros(1,N+1);
        
        for jj=1:N+1
            
            if jj~=posj
                
                pj=xmode;
                
                pj(var_id)=pp.x(jj);
                
                [pp.log_post(jj),pp.log_lik(jj)]=log_posterior_kernel(obj,pj);
                
            else
                
                pp.log_post(jj)=pp.log_post_mode;
                
                pp.log_lik(jj)=pp.log_lik_mode;
                
            end
            
        end
        
        forbidden=pp.log_post<=-obj.options.estim_penalty;
        
        pp.log_post(forbidden)=nan;
        
        pp.log_lik(forbidden)=nan;
        
        % make both curves have the same chosen type
        %------------------------------------------------
        pp.log_lik=utils.miscellaneous.apply_property(type,pp.log_post,pp.log_lik);
        
    end

end

function [tex_name,legend_]=do_one_plot(xname,db)

[~,legend_,tex_name]=utils.plot.curvature(xname,db);

end