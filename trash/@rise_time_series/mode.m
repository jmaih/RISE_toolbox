function varagout=mode(this,varargin)
this=double(this);
if size(this,3)>1
    error([mfilename,':: this operation is only defined for databases with one page'])
end
nanrows=any(isnan(this),2);
dd=this(~nanrows,:);
if isempty(dd)
    error([mfilename,':: no valid observations to compute the mode '])
end
varagout=mode(dd,varargin{:});
end
