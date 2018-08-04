function m=max(varargin)
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

if nargin==1||nargin==3
    this=varargin{1}.data;
    dim=1;
    if nargin==3
        dim=varargin{3};
    end
    m=utils.stat.max(this,[],dim);
elseif nargin==2
    if isa(varargin{1},'ts')
        this1=varargin{1}.data;
    else
        this1=varargin{1};
    end
    if isa(varargin{2},'ts')
        this2=varargin{2}.data;
    else
        this2=varargin{2};
    end
    m=utils.stat.max(this1,this2);
else
    error('number of arguments can only be 1 2 or 3')
end

end
