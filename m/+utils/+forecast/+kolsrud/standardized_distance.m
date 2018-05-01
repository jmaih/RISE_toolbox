function z=standardized_distance(y)
% standardized_distance - computes abs(y-mean(y))/std(y)
%
% ::
%
%
%   z=standardized_distance(y)
%
% Args:
%
%    - **y** [numeric] : N x T x G array, with
%      - **N** [numeric] : number of simulations/replications
%      - **T** [numeric] : sample length (time series dimension)
%      - **G** [numeric] : number of variables
%
% Returns:
%    :
%
%    - **z** [numeric] : N x T x G array of standardized distances to the
%    hypercenter
%
% Note:
%
% Example:
%
%    See also: chebyshev_distance, multivariate_chebyshev_box

% References:
% Dag Kolsrud (2015): "A time-simultaneous prediction box for a
% multivariate time series", Journal of Forecasting

s=std(y,1,1);
ybar=mean(y,1);
z=bsxfun(@rdivide,abs(bsxfun(@minus,y,ybar)),s);
end