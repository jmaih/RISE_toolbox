function this=interpolate(this,method)
if nargin<2
    method='spline';
end
tmp=double(this);
[n,p]=size(tmp);
if p>1
    error([mfilename,':: routine written to handle time series with only one variable.'])
end
w=(1:n)';
nanrows=isnan(tmp);
ti=w(~nanrows);
yi=tmp(~nanrows);
nx=numel(ti);
wnan=w(nanrows);
switch lower(method)
    case 'lagrange'
        fnan=0;
        for k=1:nx
            p=1;
            for j=[1:k-1,k+1:nx]
                p=p.*(wnan-ti(j))/(ti(k)-ti(j));
            end
            fnan=fnan+p*yi(k);
        end
    case 'spline'
        h=nan(nx-1,1);
        b=nan(nx-1,1);
        u=nan(nx-1,1);
        v=nan(nx-1,1);
        z=nan(nx,1);
        for ii=1:nx-1
            h(ii)=ti(ii+1)-ti(ii);
            b(ii)=6*(yi(ii+1)-yi(ii))/h(ii);
        end
        u(1)=2*(h(1)+h(2));
        v(1)=b(2)-b(1);
        for ii=2:nx-1
            u(ii)=2*(h(ii)+h(ii-1))-h(ii-1)^2/u(ii-1);
            v(ii)=b(ii)-b(ii-1)-h(ii-1)*v(ii-1)/u(ii-1);
        end
        z(nx)=0;
        for ii=nx-1:-1:2
            z(ii)=(v(ii)-h(ii)*z(ii+1))/u(ii);
        end
        z(1)=0;
        fnan=nan(size(wnan));
        for ii=1:nx-1
            inan=wnan>ti(ii) & wnan<ti(ii+1);
            fnan(inan)=z(ii)/(6*h(ii))*(ti(ii+1)-wnan(inan)).^3+...
                z(ii+1)/(6*h(ii))*(wnan(inan)-ti(ii)).^3+...
                (yi(ii+1)/h(ii)-z(ii+1)*h(ii)/6)*(wnan(inan)-ti(ii))+...
                (yi(ii)/h(ii)-z(ii)*h(ii)/6)*(ti(ii+1)-wnan(inan));
        end
        tmp(nanrows)=fnan;
        this=rise_time_series(this.start,tmp);
    otherwise
        error([mfilename,':: other interpolation methods not implemented yet'])
end
end
