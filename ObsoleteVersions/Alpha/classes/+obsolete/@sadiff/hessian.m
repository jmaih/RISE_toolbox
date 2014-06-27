function [Htrees,jaco_tree,code,code_expanded]=hessian(objectives,varlist,wrt)
if nargin<2
    error([mfilename,':: at least 2 input arguments should be provided'])
end
if nargin<3
    wrt=varlist;
end

% build wrt right here so that things don't get confused with the ordering
% later on

if ~isa(wrt,'sadiff')
    if ischar(wrt)
        wrt=cellstr(wrt);
    end
    prototype=sadiff('xxx');
    for irt=1:numel(wrt)
        func=wrt{irt};
        wrt{irt}=prototype;
        wrt{irt}.func=func;
        wrt{irt}.order=irt;
    end
end

[jaco_tree,code,code_expanded]=sadiff.jacobian(objectives,varlist,wrt);

Htrees=sadiff.jacobian(jaco_tree,varlist,wrt);

