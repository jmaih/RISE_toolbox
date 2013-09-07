function obj=tan(a)

if isnumeric(a.func)
    obj=planar(tan(a.func));
else
    obj=planar.multinary_operation('tan',a);
end