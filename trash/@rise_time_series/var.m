function m=var(this,varargin)
this=double(this);
n=size(this,2);
if size(this,3)>1
    error([mfilename,':: this operation is only defined for databases with one page'])
end
m=nan(1,n);
for ii=1:n
    dd=this(:,ii);
    dd=dd(~isnan(dd));
    if isempty(dd)
        error([mfilename,':: no valid observations to compute the standard deviation in column ',int2str(ii)])
    end
    m(ii)=var(dd,varargin{:});
end
end
