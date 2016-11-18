function prior=prior_setting_engine(prior,parray,id,prior_trunc,chkBounds)

if nargin<5
    
    chkBounds=[];
    
    if nargin < 4
        
        prior_trunc=[];
        
    end
    
end

if isempty(chkBounds),chkBounds=true; end

if isempty(prior_trunc),prior_trunc=1e-10; end

% for truncation
invgamma_upper_bound_truncation=10;

parray.id=id;

parray.prior_distrib=strrep(parray.prior_distrib,'_pdf','');

parray.prior_trunc=prior_trunc;

mean_stdev_flag=isnan(parray.prior_prob);

if mean_stdev_flag
    
    lqtl_mean=parray.prior_mean;
    
    uqtl_std=parray.prior_stdev;
    
else
    
    if parray.prior_prob<=0||parray.prior_prob>1
        
        error([mfilename,':: probability for parameter ',parray.name,' should be in (0,1]'])
        
    end
    
    lqtl_mean=parray.lower_quantile;
    
    uqtl_std=parray.upper_quantile;
    
end

% find the hyperparameters
[parray.a,parray.b,moments,ffinal]=distributions.(parray.prior_distrib)(lqtl_mean,uqtl_std,parray.prior_prob);

parray.prior_mean=moments.mean;

parray.prior_stdev=moments.sd;

disp([' parameter: ',upper(parray.name),', density:',upper(parray.prior_distrib),...
    ', hyperparameters: [',num2str(parray.a),' ',num2str(parray.b),'],',...
    'convergence ',num2str(ffinal)])

% get the functions of the distribution
[~,~,icdfn]=distributions.(parray.prior_distrib)();

bounds=[icdfn(prior_trunc,parray.a,parray.b),icdfn(1-prior_trunc,parray.a,parray.b)];

if isempty(parray.prior_mean)||isempty(parray.prior_stdev)
    
    disp([mfilename,'(GENTLE WARNING):: for these hyperparameters, the distribution ',...
        'does not have well-defined moments'])
    
end

the_message='';

if ismember(parray.prior_distrib,{'inv_gamma'})
    
    if bounds(2)>invgamma_upper_bound_truncation
        
        the_message=[mfilename,'(GENTLE WARNING):: upper bound of inverse gamma distribution ',...
            'truncated at ',num2str(invgamma_upper_bound_truncation)];
        
    end
    
    bounds(2) = min(bounds(2),invgamma_upper_bound_truncation);
    
end

if isfinite(parray.lower_bound),bounds(1)=max(bounds(1),parray.lower_bound);end

if isfinite(parray.upper_bound),bounds(2)=min(bounds(2),parray.upper_bound);end

% if the distribution has been truncated, say it here.
disp(the_message)

% check that the starting value is not outside the bounds
%--------------------------------------------------------
if chkBounds && (any(parray.start<bounds(1))||any(parray.start>bounds(2)))
    
    error([mfilename,':: parameter ',parray.name,' (',num2str(parray.start),') outside its bounds [',num2str(bounds),']'])
    
end

parray.lower_bound = bounds(1);

parray.upper_bound = bounds(2);

if isempty(prior)
    
    prior=parray;
    
else
    
    prior(end+1)=parray;
    
end

end