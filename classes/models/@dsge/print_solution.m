function outcell=print_solution(obj,varlist,orders,compact_form,chop)
% Print the solution of a model or vector of models
%
% ::
%
%   print_solution(obj)
%   print_solution(obj,varlist)
%   print_solution(obj,varlist,orders)
%   print_solution(obj,varlist,orders,compact_form)
%   print_solution(obj,varlist,orders,compact_form,chop)
%
% Args:
%
%    obj (rise | dsge): model object or vector of model objects
%
%    varlist (char | cellstr | {[]}): list of variables of interest
%
%    orders (numeric) : orders for which we want to
%      see the solution (default: [1:solve_order])
%
%    compact_form ({true} | false): if true, only the solution of unique
%      tuples (i,j,k) such that i<=j<=k is presented. If false, the solution
%      of all combinations is presented. i.e.
%      (i,j,k)(i,k,j)(j,i,k)(j,k,i)(k,i,j)(k,j,i)
%
%    chop ({1e-9} | numeric >0 && <1e-4) : any number, which in
%      absolute value is less than chop is set to 0.
%
% Returns:
%    :
%
%    - **outcell** [cellstr] : dummy output not important, kept for legacy
%      compatibility. Might be removed in future versions
%
% Note:
%
%    If a model is solved, say, up to 3rd order, one may still want to see the
%    first-order solution or the solution up to second-order only or any
%    combination of orders.
%

outcell_=cell(0,4);

if nargout
    
    outcell=outcell_;
    
end

if isempty(obj)
    
    return
    
end

if nargin<5
    
    chop=[];
    
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

string='';

for iobj=1:nobj
    
    if nobj>1
        
        string=int2str(iobj);
        
    end
    
    if ~isfield(obj(iobj).solution,'Tz')||isempty(obj(iobj).solution.Tz{1})
        
        fprintf(1,'\n%s\n',['MODEL ',string,...
            ' HAS NOT BEEN SOLVED AT LEAST UP TO THE FIRST ORDER']);
        
    else
        
        nsols=obj(iobj).nsols;
        
        for isol=1:nsols
            
            fprintf(1,'\nMODEL %s SOLUTION #%0.0f/%0.0f\n',string,isol,nsols);
            
            print_solution_engine(obj(iobj),varlist,compact_form,orders,chop,isol);
            
        end
        
    end
    
end

end

function print_solution_engine(obj,varlist,compact_form,orders,chop,isol)

default_vars=obj.endogenous.name(obj.endogenous.is_original & ~obj.endogenous.is_auxiliary);

defaults={
    'varlist',default_vars,@(x)ischar(x)||iscellstr(x)
    'compact_form',true,@(x)islogical(x)
    'orders',1:obj.options.solve_order,@(x)all(ismember(x,1:obj.options.solve_order))
    'chop',1e-9,@(x)isnumeric(x) && isscalar(x) && x>0 && x<1e-4}; %#ok<*ISCLSTR>

[varlist,compact_form,orders,chop]=...
    parse_arguments(defaults,...
    'varlist',varlist,'compact_form',compact_form,'orders',orders,...
    'chop',chop);

state_names={};

kept_states=[];

endo_names=obj.endogenous.name;

% get the location of the variables: can be model specific
ids=locate_variables(varlist,obj.endogenous.name);

var_names=endo_names(ids);

the_regimes=generic.describe_regimes(obj.markov_chains);

solver=parser.format_solver_name(obj.options.solver);

fprintf(1,'SOLVER :: %s\n',solver);

h=obj.markov_chains.regimes_number;

for ii=1:h
    
    myprologue={sprintf('\n%s %4.0f : %s','Regime',ii,the_regimes{ii})};
    
    build_printing_array(ii,myprologue);
    
end

    function build_printing_array(regime_index,myprologue)
        
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
        
        Tzname='T';
        
        ifact=1./cumprod(1:orders(end));
        
        for io=1:orders(end)
            
            Tzname=[Tzname,'z'];
            
            if any(io==orders)
                
                Tz=obj.solution.(Tzname){1,regime_index,isol}(ids,:).';
                
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
        
        table_displayer(bigtime,var_names,these_states,myprologue)
        
    end

end
