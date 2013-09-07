function obj=asinh(a)

if isnumeric(a.func)
    obj=planar(asinh(a.func));
else
    obj=planar.multinary_operation('asinh',a);
end
