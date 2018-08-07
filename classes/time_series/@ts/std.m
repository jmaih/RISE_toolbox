function m=std(this,varargin)
% Computes standard deviation for time series: It is an interface to MATLAB
% implementation of std, but adjusted to handle nan properly
%
% ::
%
%    varargout = std(db,varargin);
%
% Args:
%
%    db (ts object): times series object
%
%    varargin: varargin for std function in MATLAB
%
% Returns:
%    :
%
%    - **varargout**: output from std function
%

m=utils.stat.std(this.data,varargin{:});

end
