function [A,A0,SIGols,SIGml,XXi,Vols,Vml]=var_ols(y0,nlags,const)
if nargin<3
    const=true;
    if nargin<2
        nlags=4;
        if nargin<1
            error([mfilename,':: at least the data should be provided'])
        end
    end
end

if ~islogical(const)
    const=logical(const);
end

y=y0(:,nlags+1:end);
[nvar,T]=size(y);
X=nan(nvar*nlags+const,T);
if const
    X(end,:)=ones(1,T);
end

for lag=nlags:-1:1
    X((lag-1)*nvar+1:lag*nvar,:)=y0(:,(nlags+1:end)-lag);
end

A=y/X;

A0=0;
if const
    A0=A(:,end);
    A=A(:,1:end-1);
end

res=y-A*X(1:nlags*nvar,:)-A0(:,ones(1,T));

SIGols=res*res'/(T-nvar*nlags+const);
SIGml=res*res'/T;
XXi=(X*X')\eye(nvar*nlags+const);
Vols=kron(XXi,SIGols);
Vml=kron(XXi,SIGml);