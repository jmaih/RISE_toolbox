function obj=atan(a)

if isnumeric(a.func)
    obj=rise_sym(atan(a.func));
else
    obj=rise_sym.multinary_operation('atan',a);
end