function obj=uminus(a)

if ~is_zero(a)
    if ~isempty(a.args) && strcmp(a.func,'uminus')
        % Simplify -(-x) in x
        obj=rise_sym(a.args{1});
    else
        obj=rise_sym.multinary_operation('uminus',a);
    end
else
    obj=rise_sym(0);
end
