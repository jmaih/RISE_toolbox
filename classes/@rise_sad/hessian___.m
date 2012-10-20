function [H,Jacob]=hessian(string,wrt)

string=analytical_symbolic_form(string,valid_varnames,'symbolic');

Jacob=rise_sad.jacobian(string,wrt);
Jacob=Jacob(:);

n=numel(wrt);
H=cell(sum(1:n),1);
iter=0;
for irow=1:n
    Hi=rise_sad.jacobian(Jacob{irow},wrt(irow:end));
    H(iter+(1:n-irow+1))=Hi;
    iter=iter+(n-irow+1);
end

