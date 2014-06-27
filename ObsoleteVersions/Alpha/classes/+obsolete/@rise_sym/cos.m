function obj=cos(a)

if isnumeric(a.func)
    obj=rise_sym(cos(a.func));
else
    obj=rise_sym.multinary_operation('cos',a);
end