function outcell=print_solution(obj,varlist,orders,compact_form,precision,equation_format,file2save2)
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

if nargin<7
    file2save2=[];
    if nargin<6
        equation_format=[];
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
    end
end
defaults={
    'file2save2',[],@(x)ischar(x)
    };
[file2save2]=parse_arguments(defaults,'file2save2',file2save2);

nobj=numel(obj);
outcell0=cell(0,1);
string='';
for iobj=1:nobj
    if nobj>1
        string=int2str(iobj);
    end
    outcell0=[outcell0;{sprintf('\n%s',['MODEL ',string,' SOLUTION'])}]; %#ok<*AGROW>
    outcell1=print_solution_engine(obj(iobj),varlist,compact_form,precision,equation_format,orders);
    outcell0=[outcell0;outcell1];
end
if nargout
    outcell=outcell0;
else
    write_solution_out(outcell0,file2save2)
end

end

function mycell=print_solution_engine(obj,varlist,compact_form,precision,equation_format,orders)
default_vars=obj.endogenous.name(obj.endogenous.is_original & ~obj.endogenous.is_auxiliary);
defaults={
    'varlist',default_vars,@(x)ischar(x)||iscellstr(x)
    'compact_form',true,@(x)islogical(x)
    'precision',[],@(x)ischar(x)
    'equation_format',false,@(x)islogical(x)
    'orders',1:obj.options.solve_order,@(x)all(ismember(x,1:obj.options.solve_order))
    };
[varlist,compact_form,precision,equation_format,orders]=...
    parse_arguments(defaults,...
    'varlist',varlist,'compact_form',compact_form,'precision',precision,...
    'equation_format',equation_format,'orders',orders);

mycell=cell(0,1);
state_names={};
kept_states=[];
endo_names=obj.endogenous.name;
% get the location of the variables: can be model specific
ids=locate_variables(varlist,obj.endogenous.name);
var_names=endo_names(ids);

if ~isfield(obj.solution,'Tz')||isempty(obj.solution.Tz{1})
    mycell=[mycell;{sprintf('\n%s','MODEL HAS NOT BEEN SOLVED AT LEAST UP TO THE FIRST ORDER')}];
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
    solver=obj.options.solver;
    if isnumeric(solver)
        solver=int2str(solver);
    end
    mycell=[mycell;{sprintf('%s :: %s','SOLVER',solver)}];
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
        mycell=[mycell;{sprintf('\n%s %4.0f : %s','Regime',ii,the_regimes{ii})}];
        if equation_format
            B0=B{1,1};
            for icols=2:number_of_headers
                Bi=B{1,icols};
                tmp=[Bi(1,:),'='];
                for irows=2:size(Bi,1)
                    tmp=[tmp,Bi(irows,:),B0(irows,:),'+'];
                end
                tmp(isspace(tmp))=[];
                tmp=strrep(tmp,'+-','-');
                tmp=strrep(tmp,'-+','-');
                tmp=strrep(tmp,'steadystate','');
                tmp=tmp(1:end-1);
                mycell=[mycell;{sprintf('%s',tmp)}];
            end
        else
            for rr=1:nrows
                data_ii=cell(1,number_of_headers);
                for jj=1:number_of_headers
                    data_ii{jj}=B{1,jj}(rr,:);
                end
                mycell=[mycell;{sprintf(body_format,data_ii{:})}];
            end
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

function write_solution_out(mycell,file2save2)
if ~isempty(file2save2)
    fid=fopen(file2save2,'w');
else
    fid=1;
end
for irow=1:numel(mycell)
    fprintf(fid,'%s \n',mycell{irow});
end
if ~isempty(file2save2)
    fclose(fid);
end
end
