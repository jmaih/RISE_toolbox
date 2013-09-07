function [H,Jacob,auxiliary]=hessian(string,varnames,wrt)
if nargin<2
    error([mfilename,':: at least 2 input arguments should be provided'])
end
if nargin<3
    wrt=varnames;
end

Jacob=sad_forward.jacobian(string,varnames,wrt);

Jacob=Jacob(:);

n=numel(wrt);
H=cell(sum(1:n),1);
iter=0;
for irow=1:n
    Hi=sad_forward.jacobian(Jacob{irow},varnames,wrt(irow:end));
    H(iter+(1:n-irow+1))=Hi;
    iter=iter+(n-irow+1);
end

auxiliary={};

