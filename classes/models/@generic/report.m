function mytable=report(obj,destination_root,rep_items,varargin)
% Assigns the elements of interest to a rise_report.report object
%
%   - display the default inputs::
%
%        report(rise.empty(0))
%
%   - assign the reported elements in rep_items to destination_root::
%
%        report(obj,destination_root,rep_items)
%
%   - assign varargin to obj before doing the rest::
%
%        report(obj,destination_root,rep_items,varargin)
%
% Args:
%
%    obj (rise | dsge):
%    destination_root (rise_report.report): handle for the actual report
%    rep_items (char | cellstr): list of desired items to report. This list
%       can only include
%
%       - endogenous
%       - exogenous
%       - observables
%       - parameters
%       - solution
%       - estimation
%       - estimation_statistics
%       - equations
%       - code
%
% Returns:
%    :
%
%    - none (Update the report in place)
%

if isempty(obj)
    
    mydefaults=the_defaults();
    
    if nargout
        
        mytable=mydefaults;
        
    else
        
        disp_defaults(mydefaults);
        
    end
    
    return
    
end

if nargout
    
    mytable=[];
    
end

obj=set(obj,varargin{:});

if ~isa(destination_root,'rise_report.report')
    
    error('second argument should be a rise_report.report object')
    
end

nobj=numel(obj);

varlist=obj(1).options.rep_var_list;

precision=obj(1).options.rep_precision;

if isempty(rep_items)
    
    error('no item to report')
    
end

if ischar(rep_items)
    
    rep_items=cellstr(rep_items);
    
end

if ~iscellstr(rep_items)
    
    error('third argument must be a char or a cellstr')
    
end

rep_list ={'endogenous','exogenous','observables','parameters','solution',...
    'estimation','estimation_statistics','equations','code'};

bad=cellfun(@(x)~any(strcmp(x,rep_list)),rep_items,'uniformOutput',false);

if iscell(bad)
    
    bad=[bad{:}];
    
end

if any(bad)
    
    disp(rep_items(bad))
    
    error(['the items above are not valid reporting items for ',...
        'rise/dsge/svar/rfvar/stochvol objects'])
    
end

for it=1:numel(rep_items)
    
    report_engine(rep_items{it});
    
    destination_root.pagebreak()
    
end

    function report_engine(type)
        
        switch type
            
            case {'endogenous'}
                
                if nobj>1
                    
                    warning('reporting the endogenous for the first model only')
                    
                end
                
                governing=obj(1).endogenous.is_original & ...
                    ~obj(1).endogenous.is_auxiliary;
                
                this_table=[{'Model code','Description'};
                    obj(1).endogenous.name(governing)',...
                    obj(1).endogenous.tex_name(governing)'];
                
                destination_root.table('Caption','Endogenous Variables',...
                    'log',this_table,'longtable',true)
                
            case {'exogenous'}
                
                if nobj>1
                    
                    warning('reporting the exogenous for the first model only')
                    
                end
                
                this_table=[{'Model code','Description'}
                    obj(1).exogenous.name',obj(1).exogenous.tex_name'];
                
                destination_root.table('Caption','Exogenous Variables',...
                    'log',this_table,'longtable',true)
                
            case {'observables'}
                
                if nobj>1
                    
                    warning('reporting the observables for the first model only')
                    
                end
                
                this_table=[{'Model code','Description'}
                    obj(1).observables.name',obj(1).observables.tex_name'];
                
                destination_root.table('Caption','Observed Variables',...
                    'log',this_table,'longtable',true)
                
            case 'parameters'
                
                if nobj>1
                    
                    warning('reporting the parameters for the first model only')
                    
                end
                
                texnames=obj(1).parameters.tex_name';
                
                for iname=1:numel(texnames)
                    
                    underscores=texnames{iname}=='_';
                    
                    if sum(underscores)>1
                        
                        texnames{iname}=regexprep(texnames{iname},'(?<!\\)_','\\_');
                        
                    end
                    
                end
                
                this_table=[{'Model code','Description'}
                    obj(1).parameters.name',texnames];
                
                destination_root.table('Caption','Model Parameters',...
                    'log',this_table,'longtable',true)
                
            case 'solution'
                
                if isempty(varlist) && obj(1).endogenous.number>5
                    
                    warning(['Reporting the solution for perhaps too many variables. ',...
                        'All the columns might not fit a the report'])
                    
                end
                
                solution=print_solution_legacy(obj,varlist,[],[],precision);
                
                destination_root.table('Caption','Model Solution','log',solution,...
                    'longtable',true)
                
            case 'estimation'
                
                destination_root.table('Caption','Estimation Results',...
                    'log',model_estimation_results(obj),'longtable',true)
                
            case 'estimation_statistics'
                
                destination_root.table('Caption','Estimation Statistics',...
                    'log',estimation_statistics(obj))
                
            case 'equations'
                
                if nobj>1
                    
                    warning('reporting the parameters for the first model only')
                    
                end
                
                destination_root.verbatim(model_equations(obj))
                
            case 'code'
                
                model_syntax=true;
                
                model_line_numbers=true;
                
                tex_code = latex_model_file(obj,model_syntax,model_line_numbers);
                
                destination_root.text(tex_code)
                
            otherwise
                
                error(['unknown flag :: ',type])
                
        end
        
    end

end

function eqtns=model_equations(model_objects)

eqtns=cell(0,1);

dynamic=model_objects(1).equations.dynamic;

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

function model_names=get_model_names(obj)

model_names={obj.filename};

model_legends={obj.legend};

for ii=1:numel(model_legends)
    
    if ~isempty(model_legends{ii})
        
        model_names{ii}=model_legends{ii};
        
    end
    
end

end

function ptex=pick_priors_names(tex_names,model_names,threshold)

if nargin<3
    
    threshold=15;
    
end

% pick the math form first
ptex=regexprep(tex_names,'[^#]*#','');

for iname=1:numel(ptex)
    
    if length(ptex{iname})>threshold
        
        ptex{iname}=model_names{iname};
        
    end
    
    if ~any(ptex{iname}=='$')
        
        ptex{iname}=['$',ptex{iname},'$'];
        
    end
    
end

end

function estim=model_estimation_results(obj)

threshold=20;

ncases=numel(obj);

parnames=pick_priors_names({obj(1).estimation.priors.tex_name},...
    {obj(1).estimation.priors.name},...
    threshold);

ordered_names=sort(parnames);

PALL=cell(1,ncases);

for ic=1:ncases
    
    newnames=pick_priors_names({obj(ic).estimation.priors.tex_name},...
        {obj(ic).estimation.priors.name},...
        threshold);
    
    mode=num2cell(obj(ic).estimation.posterior_maximization.mode);
    
    mode_std=num2cell(obj(ic).estimation.posterior_maximization.mode_stdev);
    
    prior_prob=num2cell(vertcat(obj(ic).estimation.priors.prior_prob));
    
    plb=num2cell(vertcat(obj(ic).estimation.priors.lower_quantile));
    
    pub=num2cell(vertcat(obj(ic).estimation.priors.upper_quantile));
    
    prior_distrib={obj(ic).estimation.priors.prior_distrib};
    
    PALL{ic}=[newnames',prior_distrib',prior_prob,plb,pub,mode,mode_std];
    
    if ~isequal(parnames,newnames)
        
        ordered_names=union(ordered_names,newnames);
        
    end
    % find the bounds if they are nan...
    for irow=1:size(PALL{ic},1)
        
        if isnan(PALL{ic}{irow,3})%||PALL{ic}{irow,3}~=0.9
            
            PALL{ic}{irow,3}=0.9;
            
            bounds=distributions.find_bounds(PALL{ic}{irow,2},...
                obj(ic).estimation.priors(irow).prior_mean,...
                obj(ic).estimation.priors(irow).prior_stdev,...
                PALL{ic}{irow,3});
            
            PALL{ic}{irow,4}=bounds(1);
            
            PALL{ic}{irow,5}=bounds(2);
            
        end
        
    end
    
end

nparams=numel(ordered_names);

estim=cell(nparams+1,5+ncases);

model_names=get_model_names(obj);

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
                %                 estim{iparam+1,1}=['$',estim{iparam+1,1},'$'];
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
    
    estim=[estim,[{'mode\_std'};num2cell(obj.estimation.posterior_maximization.mode_stdev)]];
    
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
    'estimation sample'
    'number of observations '
    'number of parameters '
    'number of func. evals '
    'estimation algorithm '
    'solution algorithm '
    'start time:'
    'end time :'
    'total time:'
    };

model_names=get_model_names(obj);

for icu=1:numel(obj)
    
    this_ic={obj(icu).estimation.posterior_maximization.log_post,...
        obj(icu).estimation.posterior_maximization.log_lik,...
        obj(icu).estimation.posterior_maximization.log_prior,...
        obj(icu).estimation.posterior_maximization.log_endog_prior,...
        obj(icu).estimation.posterior_maximization.active_inequalities_number,...
        obj(icu).estimation.posterior_maximization.log_marginal_data_density_laplace,...
        [parser.any2str(obj(icu).options.estim_start_date),':',parser.any2str(obj(icu).options.estim_end_date)],...
        obj(icu).data.nobs,numel(obj(icu).estimation.priors),...
        obj(icu).estimation.posterior_maximization.funevals,...
        obj(icu).options.optimizer,...
        parser.format_solver_name(obj(icu).options.solver),...
        nan,nan,nan}';
    
    t2=obj(icu).estimation.posterior_maximization.estim_end_time;
    
    t1=obj(icu).estimation.posterior_maximization.estim_start_time;
    
    estimation_time=etime(t2,t1);
    
    hrs=floor(estimation_time/3600);
    
    secs=estimation_time-hrs*3600;
    
    mins=floor(secs/60);
    
    secs=secs-mins*60;
    
    this_ic(end-(2:-1:0))={datestr(t1),datestr(t2),[int2str(hrs),':',int2str(mins),':',int2str(secs)]};
    
    stats=[stats,[model_names(icu);this_ic]];
    
end

end

function d=the_defaults()

d={
    'rep_var_list','',@(x)ischar(x)||iscellstr(x),...
    'rep_var_list must be char or cellstr'
    
    'rep_precision','%8.6f',@(x)ischar(x),...
    'rep_precision must be char'
    };

end

function outcell=print_solution_legacy(obj,varlist,orders,compact_form,...
    precision)
% print_solution - print the solution of a model or vector of models
%
% ::
%
%
%   print_solution(obj)
%   print_solution(obj,varlist)
%   print_solution(obj,varlist,orders)
%   print_solution(obj,varlist,orders,compact_form)
%   print_solution(obj,varlist,orders,compact_form,precision)
%   print_solution(obj,varlist,orders,compact_form,precision,equation_format)
%   print_solution(obj,varlist,orders,compact_form,precision,equation_format,file2save2)
%   outcell=print_solution(obj,...)
%
% Args:
%
%    - **obj** [rise|dsge] : model object or vector of model objects
%
%    - **varlist** [char|cellstr|{[]}] : list of variables of interest
%
%    - **orders** [numeric|{[1:solve_order]}] : orders for which we want to
%      see the solution
%
%    - **compact_form** [{true}|false] : if true, only the solution of unique
%      tuples (i,j,k) such that i<=j<=k is presented. If false, the solution
%      of all combinations is presented. i.e.
%      (i,j,k)(i,k,j)(j,i,k)(j,k,i)(k,i,j)(k,j,i)
%
%    - **precision** [char|{'%8.6f'}] : precision of the numbers printed
%
%    - **equation_format** [true|{false}] : if true, the solution is presented
%      in the form of equations for each endogenous variable (not recommended)
%
%    - **file2save2** [char|{''}] : if not empty, the solution is written to a
%      file rather than printed on screen. For this to happen, print_solution
%      has to be called without ouput arguments
%
% Returns:
%    :
%
%    - **outcell** [cellstr] : If an output is requested, the solution is not
%      printed on screen or to a file.
%
% Note:
%
%    If a model is solved, say, up to 3rd order, one may still want to see the
%    first-order solution or the solution up to second-order only or any
%    combination of orders.
%
% Example:
%
%    See also:


if isempty(obj)
    
    outcell=cell(0,4);
    
    return
    
end

if nargin<5
    precision=[];
    if nargin<4
        compact_form=[];
        if nargin<3
            orders=[];
            if nargin<2
                varlist=[];
            end
        end
    end
end


nobj=numel(obj);

outcell0=cell(0,1);

string='';

for iobj=1:nobj
    
    if nobj>1
        
        string=int2str(iobj);
        
    end
    
    outcell1=print_solution_engine(obj(iobj),varlist,compact_form,precision,orders);
    
    if iobj==1
        
        ncols=size(outcell1,2);
        
        proto=cell(1,ncols);
        
        proto2=proto;
        
        outcell0=proto(1:0,:);
        
    end
    
    proto2{1}=['MODEL ',string,' SOLUTION'];
    
    outcell0=[outcell0;proto;proto2]; %#ok<*AGROW>
    
    outcell0=[outcell0;outcell1];
    
end

outcell=outcell0;

ping=cellfun(@(x)~isempty(x),outcell(:,1),'uniformOutput',true);

outcell(ping,1)=strrep(outcell(ping,1),'{','\{');

outcell(ping,1)=strrep(outcell(ping,1),'}','\}');

end

function mycell=print_solution_engine(obj,varlist,compact_form,precision,orders)

default_vars=obj.endogenous.name(obj.endogenous.is_original & ~obj.endogenous.is_auxiliary);

defaults={
    'varlist',default_vars,@(x)ischar(x)||iscellstr(x)
    'compact_form',true,@(x)islogical(x)
    'precision',[],@(x)ischar(x)
    'orders',1:obj.options.solve_order,@(x)all(ismember(x,1:obj.options.solve_order))
    };

[varlist,compact_form,precision,orders]=...
    parse_arguments(defaults,...
    'varlist',varlist,'compact_form',compact_form,'precision',precision,...
    'orders',orders);

protocell=cell(1,numel(varlist)+1);

mycell=protocell(1:0,:);

state_names={};

kept_states=[];

endo_names=obj.endogenous.name;

% get the location of the variables: can be model specific
ids=locate_variables(varlist,obj.endogenous.name);

var_names=endo_names(ids);

if ~isfield(obj.solution,'Tz')||isempty(obj.solution.Tz{1})
    
    newcell=protocell;
    
    newcell{1}='MODEL HAS NOT BEEN SOLVED AT LEAST UP TO THE FIRST ORDER';
    
    mycell=[mycell
        protocell
        newcell];
    
else
    
    the_regimes=obj.markov_chains.regimes;
    
    nregs=obj.markov_chains.regimes_number;
    
    nchains=obj.markov_chains.chains_number;
    
    chain_names=obj.markov_chains.chain_names;
    
    tmp=cell(1,nregs);
    
    for ireg=1:nregs
        
        for ichain=1:nchains
            
            new_item=[chain_names{ichain},' = ',sprintf('%0.0f',the_regimes{ireg+1,ichain+1})];
            
            if ichain==1
                
                tmp{ireg}=new_item;
                
            else
                
                tmp{ireg}=[tmp{ireg},' & ',new_item];
                
            end
            
        end
        
    end
    
    the_regimes=tmp;
    
    solver=parser.format_solver_name(obj.options.solver);
    
    newcell=protocell;
    
    newcell{1}=sprintf('%s :: %s','SOLVER',solver);
    
    mycell=[mycell
        newcell];
    
    h=obj.markov_chains.regimes_number;
    
    for ii=1:h
        
        data=build_printing_array(ii);
        
        B=utils.table.concatenate(data,precision);
        
        body_format='';
        
        % start from the end
        for bb=size(B,2):-1:1
            
            body_format=['%',int2str(B{2,bb}),'s ',body_format];
            
        end
        
        nrows=size(B{1,1},1);
        
        number_of_headers=size(B,2);
        
        newcell=protocell;
        
        newcell{1}=sprintf('\n%s %4.0f : %s','Regime',ii,the_regimes{ii});
        
        mycell=[mycell
            newcell];
        
        for rr=1:nrows
            
            data_ii=cell(1,number_of_headers);
            
            for jj=1:number_of_headers
                
                data_ii{jj}=B{1,jj}(rr,:);
                
            end
            
            mycell=[mycell
                data_ii];
            
        end
        
    end
    
end

    function data=build_printing_array(regime_index)
        
        bigtime=[
            obj.solution.ss{regime_index}(ids,:).'
            obj.solution.bgp{regime_index}(ids,:).'
            ];
        
        if regime_index==1
            
            orders=sort(orders);
            
            if ~isequal(orders,unique(orders))
                
                error('the orders specified are duplicated')
                
            end
            
            if orders(end)>obj.options.solve_order
                
                error('highest order requested exceeds solve_order')
                
            end
            
            regular=isequal(orders(:).',1:orders(end));
            
            [state_list,kept_states]=create_state_list(obj,orders,compact_form);
            
            if regular
                
                state_names=[state_names,'steady state','bal. growth']; %#ok<*AGROW>
                
                kept_states=[true(2,1);kept_states(:)]; % add steady state and bal. growth
                
            else
                
                kept_states=[false(2,1);kept_states(:)];
                
            end
            
            state_names=[state_names,state_list];
            
        end
        
        chop=1e-9;
        
        Tzname='T';
        
        ifact=1./cumprod(1:orders(end));
        
        for io=1:orders(end)
            
            Tzname=[Tzname,'z'];
            
            if any(io==orders)
                
                Tz=obj.solution.(Tzname){regime_index}(ids,:).';
                
                bigtime=[
                    bigtime
                    ifact(io)*Tz];
                
            end
            
        end
        
        if compact_form
            
            bigtime=bigtime(kept_states,:);
            
        end
        
        bigtime=full(bigtime);
        
        bigtime(abs(bigtime)<chop)=0;
        
        keep_rows=any(bigtime,2);
        
        bigtime=bigtime(keep_rows,:);
        
        % now process and print
        these_states=state_names(keep_rows);
        
        data=[{'Endo Name'},var_names
            these_states(:),num2cell(bigtime)];
        
        clear these_states bigtime
        
    end

end

