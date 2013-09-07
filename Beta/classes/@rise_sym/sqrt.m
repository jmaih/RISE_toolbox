function obj=sqrt(a)

if isnumeric(a.func)
    obj=rise_sym(sqrt(a.func));
else
    obj=rise_sym.multinary_operation('sqrt',a);
end
