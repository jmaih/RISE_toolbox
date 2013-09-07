function obj=erf(a)

if ~is_zero(a)
    obj=planar.multinary_operation('erf',a);
else
    obj=planar(0);
end
