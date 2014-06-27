function K=kurtosis(this,varargin)
this=double(this);
if size(this,3)>1
    error([mfilename,':: this operation is only defined for databases with one page'])
end
K=kurtosis(this,varargin{:});
end
