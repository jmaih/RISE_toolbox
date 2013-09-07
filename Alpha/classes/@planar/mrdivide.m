function obj=mrdivide(a,b)

if ~isa(a,'planar')
    a=planar(a);
elseif ~isa(b,'planar')
    b=planar(b);
end
if isnumeric(a.func) && isnumeric(b.func)
    if b.func==0
        error('dividing by zero not allowed')
    end
    obj=planar(mrdivide(a.func,b.func));
else
    if is_zero(a)
        obj=planar(0);
    elseif is_one(b)
        obj=a;
    elseif is_zero(b)
        error('dividing by zero not allowed')
    else
        obj=planar.multinary_operation('mrdivide',a,b);
    end
end

end