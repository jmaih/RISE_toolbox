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

eqtns=rise_sym.empty(0);
old_code=true;
if old_code
    vlist=cell2mat(strcat(varlist,','));
    vlist=vlist(1:end-1);
    
    for ifunc=1:numel(model_equations)
        myfunc=str2func(['@(',vlist,')',model_equations{ifunc},';']);
        eqtns(ifunc,1)=myfunc(args{:});
    end
else
    for ifunc=1:numel(model_equations)
        [occur,model_equations{ifunc}]=find_occurrences(model_equations{ifunc},varlist);
        % re-create the function
        var_occur=varlist(occur);
        argfun=cell2mat(strcat(var_occur,','));
        model_equations{ifunc}=str2func(['@(',argfun(1:end-1),')',model_equations{ifunc}]);
        arg_occur=args(occur);
        eqtns(ifunc,1)=model_equations{ifunc}(arg_occur{:});
    end
    
end

sizeq=size(model_equations);
% if all(sizeq>1)
    eqtns=reshape(eqtns,sizeq);
% end

end


function [occur,objectives]=find_occurrences(objectives,vlist)
if isa(objectives,'function_handle')
    objectives=func2str(objectives);
    if strcmp(objectives(1),'@')
        first_close=find(objectives==')',1,'first');
        objectives=objectives(first_close+1:end);
    end
end
if ~ischar(objectives)
    error([mfilename,':: first input must be a string or a function handle'])
end
if ischar(vlist)
    vlist=cellstr(vlist);
end
vlist=vlist(:)';
args=cell2mat(strcat(vlist,'|'));
varlist = regexp(objectives,['(?<![\w])(',args(1:end-1),')(?![\w])'],'match');
occur=ismember(vlist,varlist);
end
