function obj=plus(a,b)
% overloading of plus for coef objects
obj=coef();
obj.args={a,b};
obj.func='plus';
end
