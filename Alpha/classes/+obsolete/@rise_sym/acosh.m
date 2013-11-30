function obj=acosh(a)

if isnumeric(a.func)
    obj=rise_sym(acosh(a.func));
else
    obj=rise_sym.multinary_operation('acosh',a);
end
