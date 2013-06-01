function obj=sign(a)

if isnumeric(a.func)
    obj=rise_sym(sign(a.func));
else
    obj=rise_sym.multinary_operation('sign',a);
end

