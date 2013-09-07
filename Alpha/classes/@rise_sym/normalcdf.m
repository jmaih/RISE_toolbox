function obj=normalcdf(a,b,c)
if nargin<3
    c=rise_sym(1);
    if nargin<2
        b=rise_sym(0);
    end
end

if ~isa(a,'rise_sym')
    a=rise_sym(a);
end
if ~isa(b,'rise_sym')
    b=rise_sym(b);
end
if ~isa(c,'rise_sym')
    c=rise_sym(c);
end
if isnumeric(a.func) && isnumeric(b.func) && isnumeric(c.func)
    obj=rise_sym(normalcdf(a.func,b.func,c.func));
else
    obj=rise_sym.multinary_operation('normalcdf',a,b,c);
end

end