%  optimize_any_function: finds the optimal parameters of c logistic function
% 
%  Syntax
% 
%    gc=optimize_any_function(f,gc0,lb,ub,x1,p1,x2,p2)
% 
%  Args:
% 
%     - **f** [function handle] : function to optimize
%     - **gc0** [vector] : initial guess
%     - **lb** [vector] : lower bound of the search space
%     - **ub** [vector] : upper bound of the search space
%     - **x1,x2** [scalars] : points at which to evaluate the logistic
%     - **p1,p2** [scalars] : probability values for x1 and x2
% 
%  Returns:
% 
%     - **gc** [vector] : optimized parameters of the logistic function
%     - **fval** [scalar] : value function at gc
%     - **exitflag** [vector] : See third output of fsolve
%     - **output** [vector] : empty
% 
%  Examples:
% 
%       f=@(x,gc)1./(1+exp(-gc(1).*(x-gc(2)))) % logistic distribution
%       f=@(x,gc)1-exp(-(x/gc(1))^gc(2)) % CDF Weibull distribution
%       f=@(x,gc)1-(1-x^gc(1))^gc(2) % CDF Kumaraswamy distribution
%       f=@(x,gc)1-exp(-1/gc(1)*(x-gc(2))) % CDF exponential distribution
% 
%       gc0=zeros(2,1)
%       x1=0.17; p1=0.05; x2=0.5; p2=0.95
%       gc=optimize_any_function(f,gc0,x1,p1,x2,p2)
% 
%  See also : https://www.desmos.com/calculator
%