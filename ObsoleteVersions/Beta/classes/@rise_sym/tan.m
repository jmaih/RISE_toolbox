function obj=tan(a)

if isnumeric(a.func)
    obj=rise_sym(tan(a.func));
else
    obj=rise_sym.multinary_operation('tan',a);
end