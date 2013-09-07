function obj=mtimes(a,b)

if ~isa(a,'planar')
    a=planar(a);
elseif ~isa(b,'planar')
    b=planar(b);
end
if isnumeric(a.func) && isnumeric(b.func)
    obj=planar(mtimes(a.func,b.func));
else
    if isnumeric(a.func) && all(a.func==-1)
        obj=-b;
    elseif isnumeric(b.func) && all(b.func==-1)
        obj=-a;
    elseif is_zero(a)||is_zero(b)
        obj=planar(0);
    elseif is_one(b)
        obj=a;
    elseif is_one(a)
        obj=b;
    else
        [a,b]=planar.commute(a,b);
        obj=planar.multinary_operation('mtimes',a,b);
    end
end

end