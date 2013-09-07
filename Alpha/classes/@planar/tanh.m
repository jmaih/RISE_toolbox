function obj=acosh(a)

if isnumeric(a.func)
    obj=planar(acosh(a.func));
else
    obj=planar.multinary_operation('acosh',a);
end
