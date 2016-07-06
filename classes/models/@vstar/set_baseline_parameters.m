function obj=set_baseline_parameters(obj)

if ~isempty(obj.estimation)
    
    return
    
end

nlags=obj.nlags;

p=obj.endogenous.number;

nx=obj.deterministic.number;

obj.num_regessors=p*nlags+obj.constant+nx;

thresholds=obj.thresholds;

nthresh=numel(thresholds);

nregs=nthresh+1;

offset=0;

params=nan(10000,1);

lower_bound=params;

upper_bound=params;

pmean=params;

distribs=cell(10000,1);

distribs(1:end)={'normal'};

pstdev=ones(size(params));

is_estimated=true(size(params));

obj=load_data(obj);

B=ball_park(obj);

absB=abs(B);

low=nan; high=nan;

nx=obj.deterministic.number;

named_params=[];

for ireg_=1:nregs
    
    pp_estim=true(p,obj.num_regessors);
    
    if ireg_==1
        
        pp=B;
        %         pp=[1*eye(p),zeros(p,p*(nlags-1)+obj.constant+nx)];
        
        pnames0=cell(size(pp));
        
        for ieqtn=1:p
            
            re_offset=0;
            
            for ilag=1:nlags
                
                for ivar=1:p
                    
                    re_offset=re_offset+1;
                    
                    pnames0{ieqtn,re_offset}=sprintf('a%0.0f_%0.0f_%0.0f',ieqtn,ilag,ivar);
                    
                end
                
            end
            
            for idet=1:obj.constant+nx
                
                re_offset=re_offset+1;
                
                pnames0{ieqtn,re_offset}=sprintf('c%0.0f_%0.0f',ieqtn,idet);
                
            end
            
        end
        
        pnames=pnames0;
        
        vv=2*absB;
        
    else
        thresh=thresholds(ireg_-1);
        
        bingo=locate_variables(thresh.controlled_vars,obj.endogenous.name);
        
        pp=zeros(p,obj.num_regessors);
        
        pp_estim(:)=false;
        
        pp_estim(bingo,:)=true;
        
        pnames=strcat(pnames0,'_',int2str(ireg_-1));
        
        vv=2*absB;
        
    end
    
    add_params(1:p*obj.num_regessors,low,high,pp,pp,vv,pp_estim,'normal');
    
end

if isempty(obj.options.data)
    
    error('data are needed to initialize the parameters to estimate')
    
end

variables_locations_in_data=obj.variables_locations_in_data;

nchol=1*p*(p+1)/2;

% standard deviations
pnames=parser.create_state_list('sig',p);

stdvals=zeros(1,p);
for ieqtn=1:p
    
    pos=variables_locations_in_data.endo_id(ieqtn);
    
    stdvals(pos)=std(obj.data(pos,:));
    
end

stdvals=stdvals/5;

vv2=stdvals;

add_params(1:p,nan,nan,stdvals,stdvals,vv2,true(1,p),'inv_gamma');

% correlations
ncorr=nchol-p;

pnames=cell(1,ncorr);

iter_corr=0;

for ieqtn=1:p
    
    for ivar=1:p
        
        if ivar>=ieqtn
            
            continue
            
        end
        
        iter_corr=iter_corr+1;
        
        pnames{iter_corr}=sprintf('corr_%0.0f_%0.0f',ieqtn,ivar);
        
    end
    
end

corr_vals=zeros(1,nchol-p);

vv3=1/3;

add_params(1:nchol-p,-1,1,corr_vals,corr_vals,vv3,true(1,nchol-p),'normal');

% transition functions
%---------------------
low_g=0;

high_g=50;

g_start=0.5;

for ireg_=2:nregs
    
    thresh=thresholds(ireg_-1);
    
    pos=variables_locations_in_data.thresh_id(ireg_-1);
    
    y=obj.data(pos,:);
    
    np=thresh.np;
    
    yl=min(y); yu=max(y);
    
    cstart=thresh.threshold_priors(:,1).';
    
    cstdev=thresh.threshold_priors(:,2).';
    
    y_=linspace(yl,yu,np-1+2);
    
    y_=y_(2:end-1);
    
    bad=isnan(cstart);
    
    cstart(bad)=y_(bad);
    
    bad=isnan(cstdev);
    
    cstdev(bad)=cstart(bad);
    
    low=[low_g;yl*ones(np-1,1)];
    
    high=[high_g;yu*ones(np-1,1)];
    
    g_c=[g_start;cstart(:)];
    
    g_c_mean=g_c;
    
    g_c_mean(1)=0.5*(high_g+low_g);
    
    vv4=[sqrt(1/12*(high_g-low_g)^2);cstdev(:)];
    
    pnames=[sprintf('g_%0.0f',ireg_-1),...
        parser.create_state_list(sprintf('t%0.0f',ireg_-1),np-1)];
    
    mydistr=[{'uniform'},repmat({'normal'},1,np-1)];
    
    add_params(1:np,low,high,g_c,g_c_mean,vv4,true(1,np),mydistr);
    
end

params(offset+1:end)=[];

lower_bound(offset+1:end)=[];

upper_bound(offset+1:end)=[];

is_estimated(offset+1:end)=[];

pmean(offset+1:end)=[];

pstdev(offset+1:end)=[];

distribs(offset+1:end)=[];

tags=[];

[named_params,params,lower_bound,upper_bound,pmean,pstdev,is_estimated,distribs]=resort_all(...
    named_params,params,lower_bound,upper_bound,pmean,pstdev,is_estimated,distribs);

obj.reordering_index(tags)=1:numel(tags);

obj.parameter_values=params;

corr_par=strncmp(named_params,'corr',4);

n=numel(named_params);

obj.parameters=struct('name',{named_params},'tex_name',{named_params},...
    'number',n,...
    'is_threshold',strncmp(named_params,'t',1),...
    'is_adjust_speed',strncmp(named_params,'g',1),...
    'is_deterministic',strncmp(named_params,'c',1) & ~corr_par,...
    'is_correlation',corr_par,...
    'is_stdev',strncmp(named_params,'sig',3));

priors=struct();

for ipar=1:n
    
    if ~is_estimated(ipar)
        
        continue
        
    end
    
    if obj.options.estim_no_priors
        
        if isnan(lower_bound(ipar))
            
            lower_bound(ipar)=-2*pstdev(ipar);
            
        end
        
        if isnan(upper_bound(ipar))
            
            upper_bound(ipar)=2*pstdev(ipar);
            
        end
        
        priors.(named_params{ipar})={params(ipar),...
            max(lower_bound(ipar),-2*pstdev(ipar)),...
            min(upper_bound(ipar),2*pstdev(ipar))};
        
    else
        
        if any(isnan([lower_bound(ipar),upper_bound(ipar)]))
            
            priors.(named_params{ipar})={params(ipar),pmean(ipar),pstdev(ipar),...
                distribs{ipar}};
            
        else
            
            priors.(named_params{ipar})={params(ipar),pmean(ipar),pstdev(ipar),...
                distribs{ipar},...
                lower_bound(ipar),upper_bound(ipar)};
            
        end
        
    end
    
end

obj=setup_priors(obj,priors);

    function varargout=resort_all(varargin)
        
        varargout=varargin;
        
        for iarg=1:length(varargin)
            
            if iarg==1
                
                [varargout{iarg},tags]=sort(varargin{iarg});
                
            else
                
                varargout{iarg}=varargin{iarg}(tags);
                
            end
            
        end
        
    end

    function add_params(span,low,high,pstart,pm,pstd,pp_estim,distr)
        
        if ischar(distr)
            
            distr={distr};
            
        end
        
        stretch=offset+span;
        
        lower_bound(stretch)=low;
        
        upper_bound(stretch)=high;
        
        params(stretch)=pstart;
        
        pmean(stretch)=pm(:);
        
        pstdev(stretch)=pstd(:);
        
        is_estimated(stretch)=pp_estim(:);
        
        offset=stretch(end);
        
        named_params=[named_params,pnames(:).'];
        
        distribs(stretch)=distr;
        
    end

end
