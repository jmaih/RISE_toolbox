function obj=minus(a,b)

if ~isa(a,'rise_sym')
    a=rise_sym(a);
elseif ~isa(b,'rise_sym')
    b=rise_sym(b);
end

is_0_a=is_zero(a);
is_0_b=is_zero(b);

if is_0_b
    obj=a;
elseif is_0_a
    obj=-b;
elseif isequal(a,b)
    obj=rise_sym(0);
else
    obj=rise_sym.multinary_operation('minus',a,b);
end

end