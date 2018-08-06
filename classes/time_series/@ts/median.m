function m=median(varargin)
% Computes the median of the time series
%
% ::
%
%    m = median(db);
%    m = median(db,dim);
%
% Args:
%    db (ts object): times series object to compute the median
%    dim (integer): (optional; default 1) direction to compute the median, e.g., dim = 1 is time median of each variables, and dim = 3 would be panel median
%
% Returns:
%    :
%
%    - m (double): median values
%

this=varargin{1}.data;

dim=1;

if nargin>1
    
    dim=varargin{2};
    
end

if isvector(this) && nargin<2
    
    this=this(:);
    
end

m=utils.stat.median(this,dim);

end
