function m=mean(varargin)
% Computes the mean of the time series
%
% ::
%
%    m = mean(db);
%    m = mean(db,dim);
%
% Args:
%
%    db (ts object): times series object to compute the mean
%
%    dim (integer): (optional; default 1) direction to compute the mean,
%      e.g., dim = 1 is time mean of each variables, and dim = 3 would be
%      panel mean
%
% Returns:
%    :
%
%    - m (double): mean values
%

this=varargin{1}.data;

dim=1;

if nargin>1
    
    dim=varargin{2};
    
end

if isvector(this) && nargin<2
    
    this=this(:);
    
end

m=utils.stat.mean(this,dim);

end
