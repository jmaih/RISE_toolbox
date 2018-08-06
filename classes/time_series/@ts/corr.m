function varargout=corr(this,varargin)
% Computes correlation for time series: It is an interface to MATLAB implementation of corr
%
% ::
%
%    varargout = corr(db,varargin);
%
% Args:
%
%    db (ts object): times series object
%
%    varargin: varargin for corr function in MATLAB
%
% Returns:
%    :
%
%    varargout: output from corr function
%

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
