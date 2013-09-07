function obj=sinh(a)

if isnumeric(a.func)
    obj=planar(sinh(a.func));
else
    obj=planar.multinary_operation('sinh',a);
end