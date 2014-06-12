function out=lag_names(names,n)
if nargin<2
    n=1;
end
n=-abs(n);

out=parser.concatenate_names_number(names,n);