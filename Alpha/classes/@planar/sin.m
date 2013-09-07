function obj=sin(a)

if isnumeric(a.func)
    obj=planar(sin(a.func));
else
    obj=planar.multinary_operation('sin',a);
end