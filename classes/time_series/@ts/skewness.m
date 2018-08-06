function S=skewness(this,varargin)
% Computes skewness for time series: It is an interface to MATLAB
% implementation of skewness 
%
% ::
%
%    varargout = skewness(db,varargin);
%
% Args:
%
%    db (ts object): times series object
%
%    varargin: varargin for skewness function in MATLAB
%
% Returns:
%    :
%
%    varargout: output from skewness function
%

S=utils.stat.skewness(this.data,varargin{:});
end
