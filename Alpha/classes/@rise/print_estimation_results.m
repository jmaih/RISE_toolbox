function retcode=print_estimation_results(obj,file2save2)

if isempty(obj)
    retcode=struct();
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
    reprints=[reprints;{sprintf('\n%s\n',['MODEL ',string,' ESTIMATION RESULTS'])}];
    
    data=[{'Parameter Names',{obj(kk).estimation.priors.name}}',...
        {'distribution',upper({obj(kk).estimation.priors.prior_distrib})}',...
        {'initval',[obj(kk).estimation.priors.start]'}',...
        {'mode',obj(kk).estimation.mode}',...
        {'mode_std',obj(kk).estimation.mode_stdev}']; %	,...
%        {'mean',vertcat(obj(kk).estimated_parameters.mean)}'
    
    B=concatenate(data,precision);
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
        'log-post:',obj(kk).estimation.log_post,...
        'log-lik:',obj(kk).estimation.log_lik,...
        'log-prior:',obj(kk).estimation.log_prior,...
        'log-endog_prior',obj(kk).estimation.log_endog_prior,...
        'numberOfActiveInequalities',obj(kk).estimation.active_inequalities_number)}];
    reprints=[reprints;{sprintf('%s %8.7f %s %8.7f \n','log-MDD(Laplace)',...
        obj(kk).estimation.log_marginal_data_density.laplace,'log-MDD(MHM)',...
        obj(kk).estimation.log_marginal_data_density.modified_harmonic_mean)}];
    if ~obj.is_optimal_simple_rule_model
        reprints=[reprints;{sprintf('estimation sample is:  %s : %s (%s observations)\n',...
            parser.any2str(obj(kk).options.estim_start_date),parser.any2str(obj(kk).options.estim_end_date),...
            parser.any2str(obj(kk).data.nobs))}];
    end
    reprints=[reprints;{sprintf('solution algorithm is:  %s \n',...
        parser.any2str(obj(kk).options.solver))}];
    reprints=[reprints;{sprintf('estimation algorithm is:  %s \n',...
        parser.any2str(obj(kk).options.optimizer))}];
    reprints=[reprints;{sprintf('number of estimated parameters is:  %s \n',...
        parser.any2str(numel(obj(kk).estimation.priors)))}]; % funevals
    reprints=[reprints;{sprintf('number of function evaluations is:  %s \n',...
        parser.any2str(obj(kk).estimation.funevals))}]; % 
        t2=obj(kk).estimation.estim_end_time;
        t1=obj(kk).estimation.estim_start_time;
        estimation_time=etime(t2,t1);
        hrs=floor(estimation_time/3600);
        secs=estimation_time-hrs*3600;
        mins=floor(secs/60);
        secs=secs-mins*60;
        reprints=[reprints;{sprintf('\n %11s %s %10s %s %11s %s\n',...
            'start time:',datestr(t1),...
            'end time :',datestr(t2),...
            'total time:',[int2str(hrs),':',int2str(mins),':',int2str(secs)])}];
    reprints=[reprints;{sprintf('\n%s\n','List of issues')}];
    list_of_issues=obj(kk).list_of_issues;
    number_of_issues=size(list_of_issues,1);
    % No need to add this new issue to the list of issues coz calling this function
    % twice will add unnecessary lines to the list of issues
    if ~obj(kk).is_stable_system
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
% % if nargout
% %     if ~isempty(file2save2)
% %         error([mfilename,':: this function must be called without an extra input argument if there is an output argument'])
% %     end
% %     give_out=reprints;
% % else
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
% % end

function B=concatenate(data,precision)

B=cell(2,0);
nargs=size(data,2);

for ii=1:nargs
    if ii==1
        span=numel(data{2,1});
    end
    A=data{1,ii};
    add_on=data{2,ii};
    if iscellstr(add_on)
        add_on=char(add_on);
    end
    if isempty(add_on)
        A=char(A,repmat('--',span,1));
    elseif isnumeric(add_on)
        for jj=1:span
            A=char(A,num2str(add_on(jj),precision));
        end
    elseif ischar(A)
        A=char(A,add_on);
    else
        error([mfilename,':: unknown type'])
    end
    B=[B,{A,size(A,2)+2}']; %#ok<AGROW>
    
end