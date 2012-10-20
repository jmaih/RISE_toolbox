function J=jacobian(string,wrt)

if ischar(string)
    string=cellstr(string);
end
neq=numel(string);
if ischar(wrt)
    wrt=cellstr(wrt);
end
n=numel(wrt);
J=cell(neq,n);

for irow=1:neq
    func=string{irow};
    [J(irow,:)] = rise_sad.diff(func,wrt);
end
