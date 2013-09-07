function obj=eq(a,b)

if ~isa(a,'planar')
    a=planar(a);
elseif ~isa(b,'planar')
    b=planar(b);
end
if isnumeric(a.func) && isnumeric(b.func)
    obj=planar(eq(a.func,b.func));
else
    obj=planar.multinary_operation('eq',a,b);
end

end