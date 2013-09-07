function obj=exp(a)

if isnumeric(a.func)
    obj=rise_sym(exp(a.func));
else
    obj=rise_sym.multinary_operation('exp',a);
end