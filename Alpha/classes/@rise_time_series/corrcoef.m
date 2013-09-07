function varargout=corrcoef(this,varargin)
if ~isempty(varargin) && isa(varargin{1},'rise_time_series')
    this=this & varargin{1};
    varargin=varargin(2:end);
end
this=double(this);
if size(this,3)>1
    error([mfilename,':: this operation is only defined for databases with one page'])
end
nanrows=any(isnan(this),2);
dd=this(~nanrows,:);
if isempty(dd)
    error([mfilename,':: no valid observations to compute the correlation '])
end
[R,P,RLO,RUP]=corrcoef(dd,varargin{:});
nout=nargout;
varargout{1}=R;
if nout>1
    varargout{2}=P;
    if nout>2
        varargout{3}=RLO;
        if nout>3
            varargout{4}=RUP;
        end
    end
end
end
