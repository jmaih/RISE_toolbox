function obj=min(a,b)

if ~isa(a,'rise_sym')
    a=rise_sym(a);
elseif ~isa(b,'rise_sym')
    b=rise_sym(b);
end

obj=rise_sym.multinary_operation('min',a,b);

end