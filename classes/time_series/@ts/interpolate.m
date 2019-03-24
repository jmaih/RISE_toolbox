%--- help for ts/interpolate ---
%
%  Interpolate based on the time series to fill in the missing dates
% 
%  ::
% 
%     db = interpolate(db);
%     db = interpolate(db, method, varargin);
% 
%  Args:
% 
%     db (ts object): time series object to interpolate from
% 
%     method (string): interpolate method (default: spline). Same as the option used in MATLAB, i.e.,
% 
%         - 'nearest': nearest neighbor interpolation
%         - 'linear': linear interpolation
%         - 'spline': piecewise cubic spline interpolation (SPLINE)
%         - 'pchip': shape-preserving piecewise cubic interpolation
%         - 'cubic': same as 'pchip'
%         - 'v5cubic': the cubic interpolation from MATLAB 5, which does not
%           extrapolate and uses 'spline' if X is not equally spaced.
% 
%     varargin{1}: optional condition for extrapolation (default: no extrapolation)
% 
%        - 'extrap': extrapolate out for dates outside the dates with observations
%        - extrapval (double): use extrapval for dates outside the dates with observations
% 
%
%    Other functions named interpolate
%
%       sde/interpolate
%