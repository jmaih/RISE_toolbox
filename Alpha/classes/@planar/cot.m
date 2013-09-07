function obj=cot(a)

if isnumeric(a.func)
    obj=planar(cot(a.func));
else
    obj=planar.multinary_operation('cot',a);
end