function args=arguments(varlist,wrt)
if ischar(varlist)
    varlist=cellstr(varlist);
end

args=varlist(:)';
theargs={};
for iarg=1:numel(args)
    incidence=strcmp(args{iarg},wrt);
    if any(incidence)
        args{iarg}=rise_sym(args{iarg},theargs,incidence);
    else
        args{iarg}=rise_sym(args{iarg});
    end
end
