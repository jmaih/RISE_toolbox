function m=cov(this,varargin)
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
if size(this.data,3)>1
    error([mfilename,':: this operation is only defined for databases with one page'])
end
m=utils.stat.cov(this.data,varargin{:});
end
