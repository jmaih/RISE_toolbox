function varargout=corr(this,varargin)
if ~isempty(varargin) && isa(varargin{1},'ts')
    this=this & varargin{1};
    varargin=varargin(2:end);
end
this=this.data;
if size(this,3)>1
    error([mfilename,':: this operation is only defined for databases with one page'])
end
nout=nargout;
[varargout{1:nout}]=utils.stat.corr(this.data,varargin{:});

end
