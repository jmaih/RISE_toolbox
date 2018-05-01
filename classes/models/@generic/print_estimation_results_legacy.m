function retcode=print_estimation_results(obj,file2save2)
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

if nargin<2
    
    file2save2=[];
    
end

precision='%8.4f';

if isempty(obj(1).estimation)
    
    disp('Estimation not done!!!')
    
    return
    
end

nobj=numel(obj);

string='';

reprints=[];

for kk=1:numel(obj)
    
    if nobj>1
        
        string=int2str(kk);
        
    end
    
    is_dsge=isa(obj(kk),'dsge');
    
    reprints=[reprints;{sprintf('\n%s\n',['MODEL ',string,' ESTIMATION RESULTS'])}];
    
    data=[{'Parameter Names',{obj(kk).estimation.priors.name}}',...
        {'distribution',upper({obj(kk).estimation.priors.prior_distrib})}',...
        {'initval',[obj(kk).estimation.priors.start]'}',...
        {'mode',obj(kk).estimation.posterior_maximization.mode}',...
        {'mode_std',obj(kk).estimation.posterior_maximization.mode_stdev}']; %	,...
    %        {'mean',vertcat(obj(kk).estimated_parameters.mean)}'
    
    B=utils.table.concatenate(data,precision);
    
    body_format='\n';
    
    % start from the end
    for ii=size(B,2):-1:1
        
        body_format=['%',int2str(B{2,ii}),'s ',body_format];
        
    end
    
    nrows=size(B{1,1},1);
    
    number_of_headers=size(B,2);
    
    for ii=1:nrows
        
        data_ii=cell(1,number_of_headers);
        
        for jj=1:number_of_headers
            
            data_ii{jj}=B{1,jj}(ii,:);
            
        end
        
        reprints=[reprints;{sprintf(body_format,data_ii{:})}]; %#ok<*AGROW>
        
    end
    
    reprints=[reprints;{sprintf('\n %s %8.4f %s %8.4f %s %8.4f %s %8.4f %s %0.4g \n',...
        'log-post:',obj(kk).estimation.posterior_maximization.log_post,...
        'log-lik:',obj(kk).estimation.posterior_maximization.log_lik,...
        'log-prior:',obj(kk).estimation.posterior_maximization.log_prior,...
        'log-endog_prior',obj(kk).estimation.posterior_maximization.log_endog_prior,...
        'numberOfActiveInequalities',obj(kk).estimation.posterior_maximization.active_inequalities_number)}];
    
    reprints=[reprints;{sprintf('%s %8.7f \n','log-MDD(Laplace)',...
        obj(kk).estimation.posterior_maximization.log_marginal_data_density_laplace)}];
    
    if is_dsge && ~obj.is_optimal_simple_rule_model
        
        start_date=obj(kk).options.estim_start_date;
        
        if ~ischar(start_date)
            
            start_date=serial2date(start_date);
            
        end
        
        end_date=obj(kk).options.estim_end_date;
        
        if ~ischar(end_date)
            
            end_date=serial2date(end_date);
            
        end
        
        reprints=[reprints;{sprintf('estimation sample is:  %s : %s (%s observations)\n',...
            parser.any2str(start_date),parser.any2str(end_date),...
            parser.any2str(obj(kk).data.nobs))}];
        
    end
    
    if is_dsge
        
        reprints=[reprints;{sprintf('solution algorithm is:  %s \n',...
            parser.any2str(obj(kk).options.solver))}];
        
    end
    
    reprints=[reprints;{sprintf('estimation algorithm is:  %s \n',...
        parser.any2str(obj(kk).options.optimizer))}];
    
    reprints=[reprints;{sprintf('number of estimated parameters is:  %s \n',...
        parser.any2str(numel(obj(kk).estimation.priors)))}]; % funevals
    
    reprints=[reprints;{sprintf('number of function evaluations is:  %s \n',...
        parser.any2str(obj(kk).estimation.posterior_maximization.funevals))}]; %
    
    t2=obj(kk).estimation.posterior_maximization.estim_end_time;
    
    t1=obj(kk).estimation.posterior_maximization.estim_start_time;
    
    estimation_time=etime(t2,t1);
    
    [~,hrs,mins,secs] = utils.estim.sec2hms(estimation_time);
    
    reprints=[reprints;{sprintf('\n %11s %s %10s %s %11s %s\n',...
        'start time:',datestr(t1),...
        'end time :',datestr(t2),...
        'total time:',[int2str(hrs),':',int2str(mins),':',int2str(secs)])}];
    
    reprints=[reprints;{sprintf('\n%s\n','List of issues')}];
    
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
            
            reprints=[reprints;{sprintf('%s\n',[int2str(ii),'- ',list_of_issues{ii}])}];
            
        end
        
    else
        
        reprints=[reprints;{sprintf('%s\n','none')}];
        
    end
    
end

if ~isempty(file2save2)
    
    fid=fopen(file2save2,'w');
    
else
    
    fid=1;
    
end


for ii=1:numel(reprints)
    
    fprintf(fid,reprints{ii});
    
end

if ~isempty(file2save2)
    
    fclose(fid);
    
end

end

