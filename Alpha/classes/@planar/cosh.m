function obj=cosh(a)

if isnumeric(a.func)
    obj=planar(cosh(a.func));
else
    obj=planar.multinary_operation('cosh',a);
end