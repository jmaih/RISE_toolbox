function obj=log(a)

if isnumeric(a.func)
    if a.func<=0
        error('log is undefined for a<=0')
    end
    obj=rise_sym(log(a.func));
else
    obj=rise_sym.multinary_operation('log',a);
end