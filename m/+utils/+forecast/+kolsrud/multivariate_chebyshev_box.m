%  INTERNAL FUNCTION: constructs chebyshev boxes for multivariate-multiperiods densities
% 
%  ::
% 
%    mvcb=multivariate_chebyshev_box(y,gam,c)
% 
%  Args:
% 
%     - **y** [numeric] : N x T x G array, with
% 
%       - **N** [numeric] : number of simulations/replications
%       - **T** [numeric] : sample length (time series dimension)
%       - **G** [numeric] : number of variables
% 
%     - **gam** [scalar|vector] : percentile(s)
%     - **c** [empty|scalar|vector] : precomputed chebyshev distances
% 
%  Returns:
%     :
% 
%     - **mvcb** [2 x T x G x numel(gam)] :  array of boxes
%     - **my** [1 x T x G] :  mean across simulations
% 
%  See also:
%     - chebyshev_distance
%     - standardized_distance
% 
%