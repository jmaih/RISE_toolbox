function [lik,incr,oo_,retcode]=var_likelihood(params,y0,nlags,const)

if nargin<4
    const=true;
    if nargin<3
        nlags=4;
        if nargin<2
            error([mfilename,':: at least the params and y0 should be provided'])
        end
    end
end

n=size(y0,1);
if ~islogical(const)
    const=logical(const);
end
expect=n^2*nlags+n*const+.5*n*(n+1);
assert(numel(params)==expect,...
    ['number of parameters inconsistent with size of the VAR ',...
    'expecting ',int2str(expect),' parameters'])

lik=1e+8;
incr=[];
oo_=[];

y=y0(:,nlags+1:end);
T=size(y,2);
X=nan(n*nlags+const,T);
if const
    X(end,:)=ones(1,T);
end

for lag=nlags:-1:1
    X((lag-1)*n+1:lag*n,:)=y0(:,(nlags+1:end)-lag);
end

iter=n^2*nlags;
A=reshape(params(1:iter),n,n*nlags);
A0=0;
if const
    A0=params(iter+(1:n));
    iter=iter+n;
end

V=vartools.ivech(params(iter+(1:.5*n*(n+1))));

[flag,P,D]=ispd(V);
if ~flag
    retcode=1;
    return
end
Di=diag(1./diag(D));
% Vi=V\eye(n);
% or alternatively
% Vi=P*Di*P'; % which we don't even need
detV=prod(diag(D));
res=y-A0(:,ones(1,T))-A*X(1:n*nlags,:);

ldvp=-.5*log(detV)-.5*n*log(2*pi);

dires=diag(sqrt(diag(Di)))*P'*res;
% t=0;
% incr=nan(T,1);
% % while t<T
%     t=t+1;
%     incr(t)=dires(:,t)'*dires(:,t);
% end
incr=diag(dires'*dires);
incr=ldvp-.5*incr;
incr=-incr;
lik=sum(incr);

% % or equivalently
% lik=-(...
%     -.5*n*T*log(2*pi)...
%     -.5*T*log(det(V))...
%     -.5*trace(res'*Vi*res)...
%     );

oo_.A=A;
oo_.A0=A0;
oo_.V=V;

function [flag,uu,ss]=ispd(V)
flag=true;
% [CV,p] = chol(V);
% if p
%     flag=false;
% end
% the matrix is symmetric so I don't need the third output of the svd
[uu,ss]=svd(V);
if min(diag(ss))<1e-12
    flag=false;
end


