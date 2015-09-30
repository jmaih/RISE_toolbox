function obj=mtimes(a,b)
% overloading of mtimes for coef objects
obj=coef();
if isa(a,'coef') && isa(b,'coef')
    error('multiplication of coef x coef objects not allowed')
end
obj.args={a,b};
obj.func='mtimes';
end
