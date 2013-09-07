function [H,Jac_char,Hrs,Jac_rs]=hessian(objectives,varnames,wrt)
if nargin<2
    error([mfilename,':: at least 2 input arguments should be provided'])
end
if nargin<3
    wrt=varnames;
end

[Jac_char,Jac_rs]=sad_reverse.jacobian(objectives,varnames,wrt);
Jac_char=Jac_char(:);
Jac_rs=Jac_rs(:);

n=numel(wrt);
H=cell(sum(1:n),1);
Hrs=sad_reverse.empty(0,1);
iter=0;
for irow=1:n
    [Hi,Hrs(iter+(1:n-irow+1),1)]=sad_reverse.jacobian(Jac_rs(irow),varnames,wrt(irow:end));
    H(iter+(1:n-irow+1))=Hi;
    iter=iter+(n-irow+1);
end
