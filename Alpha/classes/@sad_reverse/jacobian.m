function [this,map,restrict_refs,mat]=jacobian(funcs,varList,wrt,vect_form)
if nargin<4
    vect_form=true;
end

if ischar(funcs)
    funcs=cellstr(funcs);
elseif iscell(funcs)
    if isa(funcs{1},'sad_reverse')
        funcs=[funcs{:}];
    end
end

if ischar(wrt)
    wrt=cellstr(wrt);
end

eqtn_nbr=numel(funcs);
wrt_nbr=numel(wrt);

varList_=sad_reverse.symbols(varList);
%%
if iscellstr(funcs)
    mytree=sad_reverse.empty(0,1);
    args=strcat(varList,',');
    args=cell2mat(args);
    for ieq=1:eqtn_nbr
        objective=['@(',args(1:end-1),')',funcs{ieq}];
        objective=str2func(objective);
        mytree(ieq,1)=objective(varList_{:});
    end
else
    mytree=funcs;
end

% vectorized form
if vect_form
    [this,index]=diff(mytree,wrt,true);
else
    % non-vectorized form
    [this,index]=diff(mytree,wrt);
end
if nargout>1
    [restrict_refs,map]=char(this,sad_reverse.create_map(index));
    if nargout>3
        mat=sad_reverse.create_matrix(map,restrict_refs,wrt_nbr,'Jac_');
    end
end
