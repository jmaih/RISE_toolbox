function var_list=initialize(var_list,wrt_list)

if ischar(var_list)
    var_list=cellstr(var_list);
end
if ischar(wrt_list)
    wrt_list=cellstr(wrt_list);
end
nwrt=numel(wrt_list);
nvlist=numel(var_list);
if numel(unique(wrt_list))~=nwrt
    error('repeated variable names in wrt_list')
end
if numel(unique(var_list))~=nvlist
    error('repeated variable names in var_list')
end

proto_incidence=false(1,nwrt);

for ivar=1:nvlist
    vname=var_list{ivar};
    loc=find(strcmp(vname,wrt_list));
    var_list{ivar}=planar(vname);
    if ~isempty(loc)
        incidence=proto_incidence;
        incidence(loc)=true;
        var_list{ivar}.incidence=sparse(incidence);
    end
end
