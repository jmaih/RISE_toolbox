function obj=log(a)

if isnumeric(a.func)
    if a.func<=0
        error('log is undefined for a<=0')
    end
    obj=planar(log(a.func));
else
    obj=planar.multinary_operation('log',a);
end