function obj=atan(a)

if isnumeric(a.func)
    obj=planar(atan(a.func));
else
    obj=planar.multinary_operation('atan',a);
end