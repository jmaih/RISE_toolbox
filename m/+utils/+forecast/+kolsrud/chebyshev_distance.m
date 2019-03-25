%  INTERNAL FUNCTION: computes chebyshev distances
% 
%  ::
% 
%    c=chebyshev_distance(y)
% 
%  Args:
% 
%     - **y** [numeric] : N x T x G array, with
% 
%       - **N** [numeric] : number of simulations/replications
%       - **T** [numeric] : sample length (time series dimension)
%       - **G** [numeric] : number of variables
% 
%  Returns:
%     :
% 
%     - **c** [N x 1 vector] : chebyshev distances
% 
%  See also:
%      - standardized_distance
%      - multivariate_chebyshev_box
% 
%