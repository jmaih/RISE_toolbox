function K=kurtosis(this,varargin)
% Computes kurtosis for time series: It is an interface to MATLAB
%   implementation of kurtosis 
%
% ::
%
%    varargout = kurtosis(db,varargin);
%
% Args:
%
%    db (ts object): times series object
%
%    varargin: varargin for kurtosis function in MATLAB
%
% Returns:
%    :
%
%    varargout: output from kurtosis function
%

K=utils.stat.kurtosis(this.data,varargin{:});
end
