function obj=acos(a)

if isnumeric(a.func)
    obj=rise_sym(acos(a.func));
else
    obj=rise_sym.multinary_operation('acos',a);
end