function eqtns=equation2rise_sym(model_equations,varlist,wrt)
if ischar(model_equations)
    model_equations=cellstr(model_equations);
end
if ischar(varlist)
    varlist=cellstr(varlist);
end
varlist=varlist(:)';

args=varlist;
theargs={};
for iarg=1:numel(args)
    incidence=strcmp(args{iarg},wrt);
    if any(incidence)
        args{iarg}=rise_sym(args{iarg},theargs,incidence);
    else
        args{iarg}=rise_sym(args{iarg});
    end
end

vlist=cell2mat(strcat(varlist,','));
vlist=vlist(1:end-1);
eqtns=rise_sym.empty(0);

for ifunc=1:numel(model_equations)
    myfunc=str2func(['@(',vlist,')',model_equations{ifunc},';']);
    eqtns(ifunc,1)=myfunc(args{:});
end
sizeq=size(model_equations);
% if all(sizeq>1)
    eqtns=reshape(eqtns,sizeq);
% end

end