function obj=cos(a)

if isnumeric(a.func)
    obj=planar(cos(a.func));
else
    obj=planar.multinary_operation('cos',a);
end