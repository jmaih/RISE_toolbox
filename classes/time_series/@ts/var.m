function m=var(this,varargin)
% Computes variance for time series: It is an interface to MATLAB implementation of var adjusted to handle nan properly
%
% ::
%
%    varargout = var(db,varargin);
%
% Args:
%    db (ts object): times series object
%    varargin: varargin for var function in MATLAB
%
% Returns:
%    :
%
%    varargout: output from var function
%

m=utils.stat.var(this.data,varargin{:});

end
