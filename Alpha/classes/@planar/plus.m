function obj=plus(a,b)

if ~isa(a,'planar')
    a=planar(a);
elseif ~isa(b,'planar')
    b=planar(b); 
end

is_0_a=is_zero(a);
is_0_b=is_zero(b);

if  ~is_0_a && ~is_0_b
    % Simplify x+(-y) in x-y
    if isequal(b.func,'uminus')
        obj=minus(a,b.args{1});
    else
        % commute to simplify the representation
        [a,b]=planar.commute(a,b);
        obj=planar.multinary_operation('plus',a,b);
    end
elseif ~is_0_a
    obj=a;
elseif ~is_0_b
    obj=b;
else
    obj=planar(0);
end

end