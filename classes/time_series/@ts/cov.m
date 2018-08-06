function m=cov(this,varargin)
% Computes covariances for time series: It is an interface to MATLAB implementation of cov adjusted to handle nan properly
%
% ::
%
%    varargout = cov(db,varargin);
%
% Args:
%    db (ts object): times series object
%    varargin: varargin for cov function in MATLAB
%
% Returns:
%    :
%
%    varargout: output from cov function
%

if ~isempty(varargin) && isa(varargin{1},'ts')
    
    this=this & varargin{1};
    
    varargin=varargin(2:end);
    
end

if size(this.data,3)>1
    
    error([mfilename,':: this operation is only defined for databases with one page'])

end

m=utils.stat.cov(this.data,varargin{:});

end
