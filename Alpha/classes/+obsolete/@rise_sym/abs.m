function obj=abs(a)

if isnumeric(a.func)
    obj=rise_sym(abs(a.func));
else
    obj=rise_sym.multinary_operation('abs',a);
end
