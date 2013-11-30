function obj=mpower(a,b)

if ~isa(a,'rise_sym')
    a=rise_sym(a);
elseif ~isa(b,'rise_sym')
    b=rise_sym(b);
end
if isnumeric(a.func) && isnumeric(b.func)
    obj=rise_sym(mpower(a.func,b.func));
else
    obj=rise_sym.multinary_operation('mpower',a,b);
end

end