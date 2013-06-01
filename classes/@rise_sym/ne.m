function obj=ne(a,b)

if ~isa(a,'rise_sym')
    a=rise_sym(a);
elseif ~isa(b,'rise_sym')
    b=rise_sym(b);
end
if isnumeric(a.func) && isnumeric(b.func)
    obj=rise_sym(ne(a.func,b.func));
else
    obj=rise_sym.multinary_operation('ne',a,b);
end

end