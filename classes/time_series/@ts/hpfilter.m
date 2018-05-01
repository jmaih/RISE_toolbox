function [hptrend,hpcycle] = hpfilter(db,lambda)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

% HP filters a collection of time series.
% 
% INPUTS 
%   db                       [double]   T*n ts object
%   lambda                   [double]   scalar, lambda parameter.
% 
% OUTPUTS 
%   hptrend                  [double]   T*n ts object, trend component of db.
%   hpcycle                  [double]   T*n ts object, cycle component of db.  
%               

y=db.data;

if nargin<2
    
    lambda = 1600;
    
end

[hptrend,hpcycle]=utils.filtering.hpfilter(y,lambda);

hptrend = ts(db.start,hptrend,db.varnames);

hpcycle = reset_data(hptrend,hpcycle);

end

