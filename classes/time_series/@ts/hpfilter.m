function [hptrend,hpcycle] = hpfilter(db,lambda)
% HP filters a collection of time series.
%
% Args:
%
%    db (ts object): time series object
%
%    lambda (double): scalar, lambda parameter. Refer to `wikipedia <https://en.wikipedia.org/wiki/Hodrick-Prescott_filter>`_
%
% Returns:
%    :
%
%       - hptrend (double): T*n ts object, trend component of db.
%       - hpcycle (double): T*n ts object, cycle component of db.
%

y=db.data;

if nargin<2
    
    lambda = 1600;
    
end

[hptrend,hpcycle]=utils.filtering.hpfilter(y,lambda);

hptrend = ts(db.start,hptrend,db.varnames);

hpcycle = set(hptrend,'data',hpcycle);

end

