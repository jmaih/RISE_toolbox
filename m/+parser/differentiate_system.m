function [derivs,numEqtns,numVars,jac_toc,original_funcs]=...
    differentiate_system(myfunc,input_list,wrt,order)

numEqtns=numel(myfunc);

myfunc=parser.remove_handles(myfunc);

myfunc=parser.symbolic_model(myfunc,input_list);

% list of symbols
symb_list=parser.collect_symbolic_list(myfunc,strcat(input_list,'_'));
% force 's0' and 's1' to enter the list
state_inputs={'s0','s1'};

input_list=input_list(:)';

for ii=1:numel(state_inputs)

    if ~any(strcmp(symb_list,state_inputs{ii}))
    
        symb_list=[symb_list,state_inputs{ii}];
    
    end
    
end
% sorting will be useful if we need to debug
symb_list=sort(symb_list);

args=splanar.initialize(symb_list,wrt);

numVars=numel(wrt);

original_funcs=myfunc;

for ifunc=1:numEqtns

    [occur,myfunc{ifunc}]=parser.find_occurrences(myfunc{ifunc},symb_list);
    % re-create the function
    var_occur=symb_list(occur);
    
    argfun=cell2mat(strcat(var_occur,','));
    
    myfunc{ifunc}=str2func(['@(',argfun(1:end-1),')',myfunc{ifunc}]);
    
    original_funcs{ifunc}=myfunc{ifunc};
    
    arg_occur=args(occur);
    
    myfunc{ifunc}=myfunc{ifunc}(arg_occur{:});

end

verbose=false;

tic

derivs=splanar.differentiate(myfunc,numVars,order,verbose);

derivs=splanar.print(derivs,input_list);

jac_toc=toc;

end
