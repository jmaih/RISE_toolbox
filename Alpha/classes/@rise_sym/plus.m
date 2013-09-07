function obj=plus(a,b)

if ~isa(a,'rise_sym')
    a=rise_sym(a);
elseif ~isa(b,'rise_sym')
    b=rise_sym(b); 
end

is_0_a=is_zero(a);
is_0_b=is_zero(b);

if  ~is_0_a && ~is_0_b
    % Simplify x+(-y) in x-y
    if isequal(b.func,'uminus')
        obj=minus(a,b.args{1});
    else
        % commute to simplify the representation
        [a,b]=rise_sym.commute(a,b);
        obj=rise_sym.multinary_operation('plus',a,b);
    end
elseif ~is_0_a
    obj=a;
elseif ~is_0_b
    obj=b;
else
    obj=rise_sym(0);
end

end