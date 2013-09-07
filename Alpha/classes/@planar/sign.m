function obj=sign(a)

if isnumeric(a.func)
    obj=planar(sign(a.func));
else
    obj=planar.multinary_operation('sign',a);
end

