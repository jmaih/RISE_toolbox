function obj=mrdivide(a,b)

if ~isa(a,'rise_sym')
    a=rise_sym(a);
elseif ~isa(b,'rise_sym')
    b=rise_sym(b);
end
if isnumeric(a.func) && isnumeric(b.func)
    if b.func==0
        error('dividing by zero not allowed')
    end
    obj=rise_sym(mrdivide(a.func,b.func));
else
    if is_zero(a)
        obj=rise_sym(0);
    elseif is_one(b)
        obj=a;
    elseif is_zero(b)
        error('dividing by zero not allowed')
    else
        obj=rise_sym.multinary_operation('mrdivide',a,b);
    end
end

end