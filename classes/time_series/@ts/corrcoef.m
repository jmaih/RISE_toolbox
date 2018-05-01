function varargout=corrcoef(this,varargin)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

if ~isempty(varargin) && isa(varargin{1},'ts')
    this=this & varargin{1};
    varargin=varargin(2:end);
end
data=this.data;
if size(data,3)>1
    error([mfilename,':: this operation is only defined for databases with one page'])
end
nout=nargout;
[varargout{1:nout}]=utils.stat.corrcoef(data,varargin{:});

% [R,P,RLO,RUP]=corrcoef(dd,varargin{:});
% nout=nargout;
% varargout{1}=R;
% if nout>1
%     varargout{2}=P;
%     if nout>2
%         varargout{3}=RLO;
%         if nout>3
%             varargout{4}=RUP;
%         end
%     end
% end
end
