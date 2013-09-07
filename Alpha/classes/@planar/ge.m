function obj=ge(a,b)

if ~isa(a,'planar')
    a=planar(a);
elseif ~isa(b,'planar')
    b=planar(b);
end
if isnumeric(a.func) && isnumeric(b.func)
    obj=planar(ge(a.func,b.func));
else
    obj=planar.multinary_operation('ge',a,b);
end

end