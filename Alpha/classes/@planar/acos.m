function obj=acos(a)

if isnumeric(a.func)
    obj=planar(acos(a.func));
else
    obj=planar.multinary_operation('acos',a);
end