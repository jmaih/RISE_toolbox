function obj=normalcdf(a,b,c)
if nargin<3
    c=planar(1);
    if nargin<2
        b=planar(0);
    end
end

if ~isa(a,'planar')
    a=planar(a);
end
if ~isa(b,'planar')
    b=planar(b);
end
if ~isa(c,'planar')
    c=planar(c);
end
if isnumeric(a.func) && isnumeric(b.func) && isnumeric(c.func)
    obj=planar(normalcdf(a.func,b.func,c.func));
else
    obj=planar.multinary_operation('normalcdf',a,b,c);
end

end