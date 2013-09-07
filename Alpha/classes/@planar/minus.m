function obj=minus(a,b)

if ~isa(a,'planar')
    a=planar(a);
elseif ~isa(b,'planar')
    b=planar(b);
end

is_0_a=is_zero(a);
is_0_b=is_zero(b);

if is_0_b
    obj=a;
elseif is_0_a
    obj=-b;
elseif isequal(a,b)
    obj=planar(0);
else
    obj=planar.multinary_operation('minus',a,b);
end

end