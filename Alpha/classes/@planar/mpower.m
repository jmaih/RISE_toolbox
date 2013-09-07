function obj=mpower(a,b)

if ~isa(a,'planar')
    a=planar(a);
elseif ~isa(b,'planar')
    b=planar(b);
end
if isnumeric(a.func) && isnumeric(b.func)
    obj=planar(mpower(a.func,b.func));
else
    obj=planar.multinary_operation('mpower',a,b);
end

end