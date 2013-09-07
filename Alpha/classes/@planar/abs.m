function obj=abs(a)

if isnumeric(a.func)
    obj=planar(abs(a.func));
else
    obj=planar.multinary_operation('abs',a);
end
