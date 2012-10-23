function [H,Jacob,myrefs]=hessian(string,varnames,wrt)
if nargin<2
    error([mfilename,':: at least 2 input arguments should be provided'])
end
if nargin<3
    wrt=varnames;
end

myrefs=[];
i_index=0;

[Jacob,myrefs,~,i_index]=rise_sad.jacobian(string,varnames,wrt,myrefs,i_index);

Jacob=Jacob(:);

n=numel(wrt);
H=cell(sum(1:n),1);
iter=0;
for irow=1:n
    [Hi,myrefs,~,i_index]=rise_sad.jacobian(Jacob{irow},varnames,wrt(irow:end),myrefs,i_index);
    H(iter+(1:n-irow+1))=Hi;
    iter=iter+(n-irow+1);
end

