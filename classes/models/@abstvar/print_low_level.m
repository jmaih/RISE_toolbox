function print_low_level(obj,caller,varlist)
% INTERNAL FUNCTION
%

colnames=stretch_variables(obj,obj.endogenous);

nendo=numel(colnames);

rownames=obj.exogenous;

if obj.constant

    rownames=['constant',rownames];

end

rownames=stretch_variables(obj,rownames);

vv=cell(numel(colnames),obj.nlags);

for ii=1:nendo

    vv(ii,:)=parser.create_state_list(colnames{ii},obj.nlags);

end

switch caller

    case 'print_solution'

        rownames=regexprep([rownames,vv(:).'],'(_)(\d+)\>','{-$2}');
        
        sList=utils.char.create_names([],'shock',nendo);

        rownames=[rownames,sList(:).'];

        stud='B';

        add_on=@(sol,id,reg)chol(sol.S(id,id,reg),'lower');

    case 'print_structural_form'

        rownames=regexprep([rownames,colnames,vv(:).'],'(_)(\d+)\>','{-$2}');

        rownames=[rownames,{'shock_stdev'}];

        stud='A';

        add_on=@(sol,id,reg)sol.S0(id,1,reg);

    otherwise

        error(['unknown caller "',caller,'"'])

end

if isempty(varlist)

    varlist=obj.endogenous;

end

% get the location of the variables: can be model specific
ids=locate_variables(varlist,obj.endogenous);

if obj.is_panel

    n=numel(ids);

    tmp=cell(n,1);

    for ii=1:n

        tmp{ii}=vec((ids(ii)-1)*obj.ng+1:ids(ii)*obj.ng);

    end

    ids=cell2mat(tmp);

end

epilogue=[];

the_regimes=generic.describe_regimes(obj.markov_chain_info);

sol=solve(obj);

for ireg=1:obj.nregs

    the_data=[sol.(stud)(ids,:,ireg),add_on(sol,ids,ireg)].';

    the_data(abs(the_data)<1e-9)=0;

    myprologue={sprintf('\n%s %4.0f : %s','Regime',ireg,the_regimes{ireg})};

    table_displayer(the_data,colnames(ids),rownames,myprologue,epilogue)

end

end