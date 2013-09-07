function eqtns=system(model_equations,varlist,wrt)
if ischar(model_equations)
    model_equations=cellstr(model_equations);
end
if ischar(varlist)
    varlist=cellstr(varlist);
end
varlist=varlist(:)';

args=rise_sym.arguments(varlist,wrt);

eqtns=rise_sym.empty(0);

for ifunc=1:numel(model_equations)
	% before creating the function, check which variables actually enter and use those as arguments
	occur=find_occurrences(model_equations{ifunc},varlist);
	if ~any(occur)
		error('no variable occurs in the function')
	end
	vlist=cell2mat(strcat(varlist(occur),','));
	vlist=vlist(1:end-1);
	
    myfunc=str2func(['@(',vlist,')',model_equations{ifunc},';']);
    eqtns(ifunc,1)=myfunc(args{:});
end
