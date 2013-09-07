function obj=uminus(a)

if ~is_zero(a)
    if ~isempty(a.args) && strcmp(a.func,'uminus')
        % Simplify -(-x) in x
        obj=planar(a.args{1});
    else
        obj=planar.multinary_operation('uminus',a);
    end
else
    obj=planar(0);
end
