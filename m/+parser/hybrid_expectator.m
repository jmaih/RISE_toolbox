function [Model_block,dic]=hybrid_expectator(Model_block,dic)

names=parser.name4hybrid_expectation();
w=names{1};
lambda=names{2};

if nargin==0
    Model_block=struct('w',w,'lambda',lambda);
    return
end

crit=ismember({w,lambda},{dic.parameters.name});
if sum(crit)==1
    error(['both parameters for hybrid expectations (',w,' and ',lambda,') have to be defined'])
end
dic.is_hybrid_expectations_model=all(crit);

if ~dic.is_hybrid_expectations_model
    return
end
const_loc=find(strcmp('const',{dic.markov_chains.name}));

aux_eqtns=cell(0,1);
ss_eqtns=cell(0,1);

for ieqtn=1:size(Model_block,1)
    eqtn=Model_block{ieqtn,1};
    new_eqtn=cell(2,0);
    for icol=1:size(eqtn,2)
        lead=eqtn{2,icol};
        if ~isempty(lead) && lead>0
            if lead>1
                error('all leads should be substituted before getting to this stage')
            end
            vname=eqtn{1,icol};
            newp=['hbe_param_',vname];
            new_var=['HBE_ENDO_',vname];
            newcol={
                '(',newp,'*', w ,'*',new_var,'+','(','1','-',newp,'*', w,')','*',vname,')'
                 [],  0 , [], 0 , [],  0    , [], [], [], [],  0 , [], 0, [], [], lead,[]
                };
            update_equation();
        else
            newcol=eqtn(:,icol);
        end
        new_eqtn=[new_eqtn,newcol]; %#ok<*AGROW>
    end
    Model_block{ieqtn,1}=new_eqtn;
end

naux=size(aux_eqtns,1);
if naux
    aux_eqtns=[repmat({nan},naux,1),aux_eqtns,repmat({'HBE eqtns'},naux,1)];
    [aux_eqtns,dic]=parser.capture_equations(dic,aux_eqtns,'model');
    Model_block=[Model_block;aux_eqtns];
    
    % auxiliary equations to be used in steady state computation
    %------------------------------------------------------------
    ss_eqtns=[repmat({nan},naux,1),ss_eqtns,repmat({'HBE eqtns'},naux,1)];
    [ss_eqtns,dic]=parser.capture_equations(dic,ss_eqtns,'steady_state_model');
    if ~isfield(dic,'auxiliary_equations')
        dic.auxiliary_equations=ss_eqtns;
    else
        dic.auxiliary_equations=[dic.auxiliary_equations;ss_eqtns];
    end
end

    function update_equation()
        if ~ismember(new_var,{dic.endogenous.name})
            added_eqtn=[new_var,'=',lambda,'*',vname,'+(1-',lambda,')*',new_var,'{-1};'];
            aux_eqtns=[aux_eqtns;added_eqtn];
            ss_eqtn=[new_var,'=',vname,';'];
            ss_eqtns=[ss_eqtns;ss_eqtn];
            dic.endogenous(end+1)=parser.listing('name',new_var);
            dic.parameters(end+1)=parser.listing('name',newp,...
                'governing_chain',const_loc);
            dic.auxiliary_variables.model{end+1}=new_var;
        end
    end

end