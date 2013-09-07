function obj=sinh(a)

if isnumeric(a.func)
    obj=rise_sym(sinh(a.func));
else
    obj=rise_sym.multinary_operation('sinh',a);
end