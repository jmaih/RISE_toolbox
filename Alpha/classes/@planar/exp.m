function obj=exp(a)

if isnumeric(a.func)
    obj=planar(exp(a.func));
else
    obj=planar.multinary_operation('exp',a);
end