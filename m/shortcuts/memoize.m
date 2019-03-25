%  memoize -- makes it possible to call a function with only the first input
% 
%  ::
% 
% 
%    f=memoize(g,[],a2,a3,...,an)
% 
%    f=memoize(g,order,a2,a3,...,an)
% 
%  Args:
% 
%     - **g** [function handle]: objective function
% 
%     - **order** [vector|empty]: helps to re-order the first input before
%     passing it to the g function
% 
%     - **ai** [anything]: additional arguments of the g function besides the
%     first one.
% 
%  Returns:
%     :
% 
%     - **f** [function handle]: memoized function such that
%     f(x)=g(x(order),varargin{:}), where "x" is the first argument of function
%     g.
% 
%  Note:
% 
%     - the number of outputs of function f is the same as that of function g.
% 
%  Example:
% 
%     See also:
%