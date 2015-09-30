function obj=mrdivide(a,b)
% overloading of mrdivide for coef objects
if ~isa(b,'double')
    error('the second argument must be a double')
end
obj=mtimes(a,1/b);
end
