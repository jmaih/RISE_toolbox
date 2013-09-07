function obj=atanh(a)

if isnumeric(a.func)
    obj=planar(atanh(a.func));
else
    obj=planar.multinary_operation('atanh',a);
end
