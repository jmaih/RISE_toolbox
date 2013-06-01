function obj=mtimes(a,b)

if ~isa(a,'rise_sym')
    a=rise_sym(a);
elseif ~isa(b,'rise_sym')
    b=rise_sym(b);
end
if isnumeric(a.func) && isnumeric(b.func)
    obj=rise_sym(mtimes(a.func,b.func));
else
    if isnumeric(a.func) && a.func==-1
        obj=-b;
    elseif isnumeric(b.func) && b.func==-1 
        obj=-a;
    elseif is_zero(a)||is_zero(b)
        obj=rise_sym(0);
    elseif is_one(b)
        obj=a;
    elseif is_one(a)
        obj=b;
    else
        [a,b]=rise_sym.commute(a,b);
        obj=rise_sym.multinary_operation('mtimes',a,b);
    end
end

end