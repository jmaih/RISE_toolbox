function obj=sqrt(a)

if isnumeric(a.func)
    obj=planar(sqrt(a.func));
else
    obj=planar.multinary_operation('sqrt',a);
end
