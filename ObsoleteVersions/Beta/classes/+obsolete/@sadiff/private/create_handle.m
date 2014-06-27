function varargout=create_handle(operCount,varargin)
varargout=varargin;
index=sprintf('%0.10g',operCount);
for ivar=1:length(varargin)
    varargout{ivar}=[varargin{ivar},'_',index,'_'];
end