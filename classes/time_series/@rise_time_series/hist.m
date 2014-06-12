function h=hist(this,varargin)
this=double(this);
if size(this,3)>1
    error([mfilename,':: this operation is only defined for databases with one page'])
end
nanrows=any(isnan(this),2);
dd=this(~nanrows,:);
if isempty(dd)
    error([mfilename,':: no valid observations to compute the mode '])
end
if isempty(varargin)
    hist(dd);
else
    nbins=varargin{1};
    if isnumeric(nbins)
    hist(dd,nbins);
    else
        hist(dd);
%         warning('second argument of hist is not numeric and was ignored')
    end
end
if nargout
    h=gca;
end
end
