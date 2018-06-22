function m=sum(this,varargin)
% Overloaded sum function for ts object
%

m=utils.stat.sum(this.data,varargin{:});

end
