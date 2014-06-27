function obj=asinh(a)

if isnumeric(a.func)
    obj=rise_sym(asinh(a.func));
else
    obj=rise_sym.multinary_operation('asinh',a);
end
