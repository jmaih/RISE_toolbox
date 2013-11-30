function obj=erf(a)

if ~is_zero(a)
    obj=rise_sym.multinary_operation('erf',a);
else
    obj=rise_sym(0);
end
