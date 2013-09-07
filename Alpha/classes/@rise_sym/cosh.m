function obj=cosh(a)

if isnumeric(a.func)
    obj=rise_sym(cosh(a.func));
else
    obj=rise_sym.multinary_operation('cosh',a);
end