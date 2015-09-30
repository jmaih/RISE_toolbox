function obj=minus(a,b)
% overloading of minus for coef objects
obj=coef();
obj.args={a,b};
obj.func='minus';
end
