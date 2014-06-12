function [links,parameter_values]=parameters_to_matrices(param_template,...
    param_names,numberOfRegimes)
if nargin<3
    numberOfRegimes=1;
end
% link the parameters to the structural matrices
%-----------------------------------------------
ncell=size(param_template,2);
links=cell(1,ncell);
for ii=1:ncell
    themat=param_template{2,ii};
    header=param_template{1,ii};
    [jj,kk]=find(isnan(themat));
    if isempty(jj)
        continue
    end
    plist=cellstr(strcat(header,'_',int2str(jj),'_',int2str(kk)));
    ploc=locate_variables(plist,param_names);
    links{ii}=sub2ind(size(themat),jj,kk)+ploc*1i;
end

nparams=numel(param_names);
parameter_values=nan(nparams,numberOfRegimes);