function retcode=print_estimation_results(obj)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:


if isempty(obj)
    
    retcode=cell(0,4);
    
    return
    
end

if nargout
    
    retcode=0;
    
end

if isempty(obj(1).estimation)
    
    disp('Estimation not done!!!')
    
    return
    
end

nobj=numel(obj);

string='';

for kk=1:numel(obj)
    
    if nobj>1
        
        string=int2str(kk);
        
    end
    
    is_dsge=isa(obj(kk),'dsge');
    
    myprologue={sprintf('\n%s',['MODEL ',string,' ESTIMATION RESULTS'])};
    
    mode_std=obj(kk).estimation.posterior_maximization.mode_stdev;
        
    the_data={
        char(upper({obj(kk).estimation.priors.prior_distrib})),...
        [obj(kk).estimation.priors.start].',...
        obj(kk).estimation.posterior_maximization.mode};
    
    colnames={'distribution','initval','mode'};
    
    is_hessian=~all(isnan(mode_std));
    
    if is_hessian
        
        the_data=[the_data,{mode_std}];
        
        colnames=[colnames,{'mode_std'}];
        
    end
    
    rownames={obj(kk).estimation.priors.name};
    
    myepilogue={
        sprintf('%s %8.4f %s %8.4f %s %8.4f %s %8.4f %s %0.4g',...
        'log-post:',obj(kk).estimation.posterior_maximization.log_post,...
        'log-lik:',obj(kk).estimation.posterior_maximization.log_lik,...
        'log-prior:',obj(kk).estimation.posterior_maximization.log_prior,...
        'log-endog_prior',obj(kk).estimation.posterior_maximization.log_endog_prior,...
        'numberOfActiveInequalities',obj(kk).estimation.posterior_maximization.active_inequalities_number)
        };
    
    if is_hessian
        
        myepilogue=[myepilogue
            {sprintf('%s %8.7f','log-MDD(Laplace)',...
            obj(kk).estimation.posterior_maximization.log_marginal_data_density_laplace)}
            ];
        
    end
    
    if is_dsge && ~obj(kk).is_optimal_simple_rule_model
        
        start_date=obj(kk).options.estim_start_date;
        
        if ~ischar(start_date)
            
            start_date=serial2date(start_date);
            
        end
        
        end_date=obj(kk).options.estim_end_date;
        
        if ~ischar(end_date)
            
            end_date=serial2date(end_date);
            
        end
        
        myepilogue=[myepilogue
            {sprintf('estimation sample is:  %s : %s (%s observations)',...
            parser.any2str(start_date),parser.any2str(end_date),...
            parser.any2str(obj(kk).data.nobs))}];
        
    end
    
    if is_dsge
        
        myepilogue=[myepilogue
            {sprintf('solution algorithm is:  %s',...
            parser.any2str(obj(kk).options.solver))}];
        
    end
    
    myepilogue=[myepilogue;{sprintf('estimation algorithm is:  %s',...
        parser.any2str(obj(kk).options.optimizer))}];
    
    myepilogue=[myepilogue;{sprintf('number of estimated parameters is:  %s',...
        parser.any2str(numel(obj(kk).estimation.priors)))}]; % funevals
    
    myepilogue=[myepilogue;{sprintf('number of function evaluations is:  %s',...
        parser.any2str(obj(kk).estimation.posterior_maximization.funevals))}]; %
    
    t2=obj(kk).estimation.posterior_maximization.estim_end_time;
    
    t1=obj(kk).estimation.posterior_maximization.estim_start_time;
    
    estimation_time=etime(t2,t1);
    
    [~,hrs,mins,secs] = utils.estim.sec2hms(estimation_time);
    
    myepilogue=[myepilogue;{sprintf('\n%11s %s %10s %s %11s %s',...
        'start time:',datestr(t1),...
        'end time :',datestr(t2),...
        'total time:',[int2str(hrs),':',int2str(mins),':',int2str(secs)])}];
    
    myepilogue=[myepilogue;{sprintf('\n%s\n','List of issues')}];
    
    list_of_issues=obj(kk).list_of_issues;
    
    number_of_issues=size(list_of_issues,1);
    
    % No need to add this new issue to the list of issues coz calling this function
    % twice will add unnecessary lines to the list of issues
    if is_dsge && ~obj(kk).is_stable_system
        
        number_of_issues=number_of_issues+1;
        
        list_of_issues{number_of_issues}='The system is unstable';
        
    end
    
    if number_of_issues
        
        for ii=1:number_of_issues
            
            myepilogue=[myepilogue;{sprintf('%s',[int2str(ii),'- ',list_of_issues{ii}])}];
            
        end
        
    else
        
        myepilogue=[myepilogue;{sprintf('%s','none')}]; %#ok<*AGROW>
        
    end
    
    table_displayer(the_data,colnames,rownames,myprologue,myepilogue);
    
end


end

