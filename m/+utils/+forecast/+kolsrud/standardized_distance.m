%  INTERNAL FUNCTION: computes abs(y-mean(y))/std(y)
% 
%  ::
% 
%    z=standardized_distance(y)
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
%     - **z** [numeric] : N x T x G array of standardized distances to the
%       hypercenter
% 
%  See also:
%     - chebyshev_distance
%     - multivariate_chebyshev_box
% 
%