function obj=asin(a)

if isnumeric(a.func)
    obj=planar(asin(a.func));
else
    obj=planar.multinary_operation('asin',a);
end