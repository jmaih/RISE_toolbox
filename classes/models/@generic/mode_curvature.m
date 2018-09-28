function db=mode_curvature(obj,varlist,N,type,dbin)
% Checks the curvature at the posterior mode
%
% ::
%
%   db = mode_curvature(obj)
%
%   db = mode_curvature(obj,varlist)
%
%   db = mode_curvature(obj,varlist,N)
%
%   db = mode_curvature(obj,varlist,N,type)
%
% Args:
%
%    obj (rise | dsge | rfvar | svar): model object
%
%    varlist (char | cellstr | empty): list of parameters for which we want
%      to check the curvature
%
%    N ({20} | integer): Number of grid points
%
%    type ({'max'} | 'min' | 'range'): normalization of the log-posterior
%      and the log-likelihood.
%
%    dbin (struct|empty): structure containing the information to plot the
%      curvature. Each field is the name of a particular parameter. This is
%      to avoid a costly recomputation of db
%
% Returns:
%    :
%
%    - **db** [struct|cell array|vector]: structure containing the  
%      information to plot the curvature. Each field is the name of a 
%      particular parameter. Alternatively, when dbin is not empty, db is a 
%      handle to the plots.
%
% Note:
%
%    - when no output is requested, plots are made but not saved.
%
%    - one way to plot the curvatures from the output is to use the function
%      utils.plot.curvature
%
% See also:
%    - utils.plot.curvature
%

if isempty(obj)
    
    if nargout>1
        
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
        
    end
    
    db=cell(0,4);
    
    return
    
end

if nargin<5
    
    dbin=[];
    
    if nargin<4
        
        type=[];
        
        if nargin<3
            
            N=[];
            
            if nargin<2
                
                varlist=[];
                
            end
            
        end
        
    end
    
end

if isempty(N)
    
    N=20;
    
end

if isempty(type)
    
    type='max';
    
end

nobj=numel(obj);

if nobj>1
    
    db_=cell(1,nobj);
    
    if isempty(dbin)
        
        dbin=db_;
        
    end
    
    nworkers=utils.parallel.get_number_of_workers();
    
    plot_style=nargout==0;
    
    if plot_style
        
        parfor(iobj=1:nobj,nworkers)
            
            mode_curvature(obj(iobj),varlist,N,type,dbin{iobj});
            
        end
        
    else
        
        parfor(iobj=1:nobj,nworkers)
            
            db_{iobj}=mode_curvature(obj(iobj),varlist,N,type,dbin{iobj});
            
        end
        
        db=db_;
        
    end
    
    return
    
end

if ~isempty(dbin)
    
    AllNames=fieldnames(dbin);
    
else
    
    AllNames={obj.estimation.priors.name};
    
end

if isempty(varlist)
    
    varlist=AllNames;
    
end

if ischar(varlist)
    
    varlist=cellstr(varlist);
    
end

list_locs=locate_variables(varlist,AllNames);

varnames=AllNames(list_locs); %<---{obj.estimation.priors(list_locs).name};

if ~isempty(dbin)
    
    db_=dbin;
    
else
    
    db_=do_it();
    
end

plotit=nargout==0||~isempty(dbin);

if plotit
    % plot the data
    %--------------
    r0=obj.options.graphics(1);
    
    c0=obj.options.graphics(2);
    
    titel='Curvature around the found mode';
    
    h=utils.plot.multiple(@(xname)do_one_plot(xname,db_),...
        varnames,titel,r0,c0,...
        'FontSize',11,'FontWeight','normal');
    
    if nargout
        
        db=h;
        
    end
    
else
    
    db=db_;
    
end

    function db_=do_it()
        
        LB=vertcat(obj.estimation.priors.lower_bound);
        
        UB=vertcat(obj.estimation.priors.upper_bound);
        
        xmode=obj.estimation.posterior_maximization.mode;
        
        obj.options.kf_filtering_level=0; % do not filter
        
        npar=numel(list_locs);
        
        db_=cell(2,npar);
        
        mainfunc=@do_one_parameter;
        
        % disable key warnings
        %------------------------
        warnstate=utils.estim.warnings_disable();
        
        nworkers=utils.parallel.get_number_of_workers();
        
        parfor(ipar=1:npar,nworkers)
            
            var_id=list_locs(ipar);
            
            vname=varnames{ipar};
            
            db_(:,ipar)={vname,mainfunc(var_id)}';
            
        end
        % re-enable warnings
        %-------------------
        utils.estim.warnings_enable(warnstate);
        
        fields=db_(1,:);
        
        db_=cell2struct(db_(2,:),fields,2);
        
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

end

function [tex_name,legend_]=do_one_plot(xname,db)

[~,legend_,tex_name]=utils.plot.curvature(xname,db);

end