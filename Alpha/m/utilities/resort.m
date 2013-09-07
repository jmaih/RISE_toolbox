function varargout=resort(varargin)
if nargout~=nargin
    error([mfilename,':: number of inputs and outputs should be the same'])
end
[varargout{1},tags]=sort(varargin{1});
for ii=2:nargout
    varargout{ii}=varargin{ii}(:,tags);
end

end
