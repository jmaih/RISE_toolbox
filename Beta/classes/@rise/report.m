function mytable=report(obj,varargin)
% rep_type ='varendo','varexo','varobs','parameters','solution'
% 'estimation','estimation_statistics','equations'

Defaults=struct('rep_type','varendo',...
    'rep_var_list','',...
    'rep_precision','%8.6f');

if isempty(obj)
    mytable=Defaults;
    return
end

obj=set_options(obj,varargin{:});

nobj=numel(obj);
type=obj(1).options.rep_type;
varlist=obj(1).options.rep_var_list;
precision=obj(1).options.rep_precision;
switch type
    case 'varendo'% endogenous
        if nobj>1
            warning('reporting the endogenous for the first model only')
        end
        mytable=[{'Model code','Description'};
            {obj(1).varendo.name}',{obj(1).varendo.tex_name}'];
    case 'varexo'
        if nobj>1
            warning('reporting the exogenous for the first model only')
        end
        mytable=[{'Model code','Description'}
            {obj(1).varexo.name}',{obj(1).varexo.tex_name}'];
    case 'varobs'
        if nobj>1
            warning('reporting the observables for the first model only')
        end
        mytable=[{'Model code','Description'}
            {obj(1).varobs.name}',{obj(1).varobs.tex_name}'];
    case 'parameters'
        if nobj>1
            warning('reporting the parameters for the first model only')
        end
        texnames={obj(1).parameters.tex_name}';
        for iname=1:numel(texnames)
            underscores=texnames{iname}=='_';
            if sum(underscores)>1
                texnames{iname}=regexprep(texnames{iname},'(?<!\\)_','\\_');
            end
        end
        mytable=[{'Model code','Description'}
            {obj(1).parameters.name}',...
            strcat('$',texnames,'$')];
    case 'solution'
        if isempty(varlist) && obj(1).NumberOfEndogenous(2)>5
            warning(['Reporting the solution for perhaps too many variables. ',...
                'All the columns might not fit a the report'])
        end
        mytable=print_solution(obj,varlist,precision);
    case 'estimation'
        mytable=model_estimation_results(obj);
    case 'estimation_statistics'
        mytable=estimation_statistics(obj);
    case 'equations'
        if nobj>1
            warning('reporting the parameters for the first model only')
        end
        mytable=model_equations(obj);
    otherwise
        error(['unknown flag :: ',type])
end

end

function eqtns=model_equations(model_objects)
eqtns=cell(0,1);
dynamic={model_objects(1).equations.dynamic};

for ieq=1:numel(dynamic)
    myeq=['EQ',int2str(ieq),': ',dynamic{ieq}];
    equality=strfind(myeq,'=');
    if isempty(equality)
        myeq=strrep(myeq,';','=0;');
    end
    if length(myeq)>89
        neweq='';
        strfill='';
        while ~isempty(myeq)
            [token,remain] = strtok(myeq,'+-/*');
            strfill=[strfill,token];
            if ~isempty(remain)
                strfill=[strfill,remain(1)]; %#ok<*AGROW>
            end
            myeq=remain(2:end);
            if length(strfill)>86
                if ~isempty(remain(2:end))
                    strfill=[strfill,'...'];
                end
                neweq=char(neweq,strfill);
                if ~isempty(remain(2:end))
                    strfill='	';
                else
                    strfill='';
                end
            end
        end
        if ~isempty(strfill)
            neweq=char(neweq,strfill);
        end
        myeq=neweq(2:end,:);
        clear neweq
    end
    myeq=cellstr(myeq);
    eqtns=[eqtns;myeq(:);{' '}];
end
end

function estim=model_estimation_results(obj)
ncases=numel(obj);
type_name='tex_name';
parnames= {obj(1).estimated_parameters.(type_name)};
ordered_names=sort(parnames);
PALL=cell(1,ncases);
for ic=1:ncases
    newnames={obj(ic).estimated_parameters.(type_name)};
    mode=num2cell(vertcat(obj(ic).estimated_parameters.mode));
    mode_std=num2cell(vertcat(obj(ic).estimated_parameters.mode_std));
    prior_prob=num2cell(vertcat(obj(ic).estimated_parameters.interval_probability));
    plb=num2cell(vertcat(obj(ic).estimated_parameters.plb));
    pub=num2cell(vertcat(obj(ic).estimated_parameters.pub));
    prior_distrib={obj(ic).estimated_parameters.distribution};
    PALL{ic}=[newnames',prior_distrib',prior_prob,plb,pub,mode,mode_std];
    if ~isequal(parnames,newnames)
        ordered_names=union(ordered_names,newnames);
    end
end
nparams=numel(ordered_names);

estim=cell(nparams+1,5+ncases);
model_names={obj.filename};
estim(1,:)=[{'parameter','Prior distr','Prior prob','low','high'},model_names];
for iparam=1:nparams
    name_in=false;
    for imod=1:ncases
        if isempty(PALL{imod})
            estim{iparam+1,5+imod}='--';
        else
            if ~name_in
                param_info=PALL{imod}(1,:);
                pname=param_info{1};
                name_in=true;
                estim(iparam+1,1:5)=param_info(1:5);
                estim{iparam+1,5+imod}=param_info{6};
                PALL{imod}(1,:)=[];
                % add the dollars
                estim{iparam+1,1}=['$',estim{iparam+1,1},'$'];
                underscores=estim{iparam+1,1}=='_';
                if sum(underscores)>1
                    estim{iparam+1,1}=regexprep(estim{iparam+1,1},'(?<!\\)_','\\_');
                end
                continue
            end
            loc=find(strcmp(pname,PALL{imod}(:,1)));
            if isempty(loc)
                estim{iparam+1,5+imod}='--';
            else
                estim{iparam+1,5+imod}=PALL{imod}{loc,6};
                PALL{imod}(loc,:)=[];
            end
        end
    end
end
if ncases==1
    estim=[estim,[{'mode\_std'};num2cell(vertcat(obj.estimated_parameters.mode_std))]];
end

end

function stats=estimation_statistics(obj)
stats={
    ''
    'log-post:'
    'log-lik:'
    'log-prior:'
    'log-endog_prior'
    'number of active inequalities'
    'log-MDD(Laplace)'
    'log-MDD(MHM)'
    'estimation sample'
    'number of observations '
    'number of parameters '
    'estimation algorithm '
    'solution algorithm '
    'start time:'
    'end time :'
    'total time:'
    };
model_names={obj.filename};
for icu=1:numel(obj)
    this_ic={obj(icu).log_post,obj(icu).log_lik,obj(icu).log_prior,...
        obj(icu).log_endog_prior,obj(icu).numberOfActiveInequalities,...
        obj(icu).log_mdd_laplace,obj(icu).log_mdd_mhm,...
        [obj(icu).options.estim_start_date,':',obj(icu).options.estim_end_date],...
        numel(obj(icu).varobs(1).value),numel(obj(icu).estimated_parameters),...
        obj(icu).options.optimizer,...
        obj(icu).options.solver,nan,nan,nan}';
    if isfield(obj(icu).options,'estim_end_time') && ~isempty(obj(icu).options.estim_end_time)
        t2=obj(icu).options.estim_end_time;
        t1=obj(icu).options.estim_start_time;
        estimation_time=etime(t2,t1);
        hrs=floor(estimation_time/3600);
        secs=estimation_time-hrs*3600;
        mins=floor(secs/60);
        secs=secs-mins*60;
        this_ic(end-(2:-1:0))={datestr(t1),datestr(t2),[int2str(hrs),':',int2str(mins),':',int2str(secs)]};
    end
    stats=[stats,[model_names(icu);this_ic]]; %#ok<AGROW>
end
end
