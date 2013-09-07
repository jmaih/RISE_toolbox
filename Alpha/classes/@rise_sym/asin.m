function obj=asin(a)

if isnumeric(a.func)
    obj=rise_sym(asin(a.func));
else
    obj=rise_sym.multinary_operation('asin',a);
end