function obj=sin(a)

if isnumeric(a.func)
    obj=rise_sym(sin(a.func));
else
    obj=rise_sym.multinary_operation('sin',a);
end