function obj=atanh(a)

if isnumeric(a.func)
    obj=rise_sym(atanh(a.func));
else
    obj=rise_sym.multinary_operation('atanh',a);
end
