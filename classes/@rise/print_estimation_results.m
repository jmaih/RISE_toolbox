function print_estimation_results(obj,file2save2)

if nargin<2
    file2save2=[];
end
precision='%8.4f';

nobj=numel(obj);
string='';
for kk=1:numel(obj)
    if nobj>1
        string=int2str(kk);
    end
    reprints={sprintf('\n%s\n',['MODEL ',string,' ESTIMATION RESULTS'])};
    
    data=[{'Parameter Names',{obj(kk).estimated_parameters.name}}',...
        {'distribution',upper({obj(kk).estimated_parameters.distribution})}',...
        {'initval',vertcat(obj(kk).estimated_parameters.startval)}',...
        {'mode',vertcat(obj(kk).estimated_parameters.mode)}',...
        {'mode_std',vertcat(obj(kk).estimated_parameters.mode_std)}',...
        {'mean',vertcat(obj(kk).estimated_parameters.mean)}'];
    
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
        'log-post:',obj(kk).log_post,...
        'log-lik:',obj(kk).log_lik,...
        'log-prior:',obj(kk).log_prior,...
        'log-endog_prior',obj(kk).log_endog_prior,...
        'numberOfActiveInequalities',obj(kk).numberOfActiveInequalities)}];
    reprints=[reprints;{sprintf('%s %8.7f %s %8.7f \n','log-MDD(Laplace)',...
        obj(kk).log_mdd_laplace,'log-MDD(MHM)',obj(kk).log_mdd_mhm)}];
    reprints=[reprints;{sprintf('estimation sample is:  %s : %s (%s observations)\n',...
        any2str(obj(kk).options.estim_start_date),any2str(obj(kk).options.estim_end_date),...
        any2str(numel(obj(kk).varobs(1).value)))}];
    reprints=[reprints;{sprintf('estimation algorithm is:  %s \n',...
        any2str(obj(kk).options.optimizer))}];
    reprints=[reprints;{sprintf('number of estimated parameters is:  %s \n',...
        any2str(numel(obj(kk).estimated_parameters)))}];
    reprints=[reprints;{sprintf('solution algorithm is:  %s \n',...
        any2str(obj(kk).options.solver))}];
    if isfield(obj(kk).options,'estim_end_time') && ~isempty(obj(kk).options.estim_end_time)
        t2=obj(kk).options.estim_end_time;
        t1=obj(kk).options.estim_start_time;
        estimation_time=etime(t2,t1);
        hrs=floor(estimation_time/3600);
        secs=estimation_time-hrs*3600;
        mins=floor(secs/60);
        secs=secs-mins*60;
        reprints=[reprints;{sprintf('\n %11s %s %10s %s %11s %s\n',...
            'start time:',datestr(t1),...
            'end time :',datestr(t2),...
            'total time:',[int2str(hrs),':',int2str(mins),':',int2str(secs)])}];
    end
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

function str=any2str(whatever)
if ischar(whatever)
    str=whatever;
elseif isnumeric(whatever)
    str=num2str(whatever);
elseif isa(whatever,'function_handle')
    str=func2str(whatever);
elseif iscell(whatever)
    str=any2str(whatever{1});
elseif isa(whatever,'rise_date')
    str=any2str(whatever.date);
else
    error([mfilename,':: cannot convert elements of the ',class(whatever),' class to a string'])
end


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