%  INTERNAL FUNCTION: Approximate empirical cumulative distribution function
% 
%  ::
% 
%     [f, xx] = empirical_cdf(x,lu,ub,N);
%     [f, xx] = empirical_cdf(x,lu,ub);
%     [f, xx] = empirical_cdf(x,lu);
%     [f, xx] = empirical_cdf(x);
% 
%  Args:
%     x (vector of double): data points to create empirical CDF of
%     lb (double): lower bound of the range (default min(x))
%     ub (double): upper bound of the range (default max(x))
%     N (integer): number of bins/grid points to use to approximate the empirical CDF (default 250)
% 
%  Returns:
%     :
%     xx (N x 1 double): knot points of the emprical CDF (uniformly distributed between lb and ub)
%     f (N x 1 double): value of the CDF at the given point
% 
%