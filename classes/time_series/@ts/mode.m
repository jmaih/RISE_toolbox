function varargout=mode(this,varargin)
% Computes mode for the time series: It is an interface to MATLAB
%   implementation of mode 
%
% ::
%
%    varargout = mode(db,varargin);
%
% Args:
%    db (ts object): times series object
%    varargin: varargin for mode function in MATLAB
%
% Returns:
%    :
%
%    varargout: output from mode function
%

nout=nargout;

[varargout{1:nout}]=utils.stat.mode(this.data,varargin{:});

end
