function varargout=corr(this,varargin)
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

bad=any(isnan(data),2);

if any(bad)
    
    data(bad,:)=[];
        
    disp([int2str(sum(bad)),' nan rows above discarded'])
    
end

[varargout{1:nout}]=utils.stat.corr(data,varargin{:});

end
