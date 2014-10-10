function [results,stats]=ols(y0,xdata,nlags,const,expand_results,order)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

if nargin<6
    order=5;
    if nargin<5
        expand_results=false;
        if nargin<4
            const=true;
            if nargin<3
                nlags=4;
                if nargin<2
                    xdata=[];
                    if nargin<1
                        error([mfilename,':: at least the data should be provided'])
                    end
                end
            end
        end
    end
end

if ~islogical(const)
    const=logical(const);
end

nx=size(xdata,1)+const;
y=y0(:,nlags+1:end);
[nvar,T]=size(y);
X=nan(nvar*nlags+nx,T);
if nx-const>0
    X(nvar*nlags+(1:nx-const),:)=xdata(:,nlags+1:end);
end
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
end

resid=y-A*X;

RSS=(resid*resid');
SIGols=RSS/(T-nvar*nlags+const);
results=struct('A',A,'const',A0,'residuals',resid,'SIGols',SIGols,'RSS',RSS,'A0',A0);
k=nvar*nlags+nx;
if expand_results || nargout>1
    results.SIGml=RSS/T; % maximum likelihood estimator
    results.XXi=(X*X')\eye(k);
    results.Vols=kron(results.XXi,SIGols);
    results.Vml=kron(results.XXi,results.SIGml);
end

if nargout>1
    
    ybar=mean(y,2);
    TSS=bsxfun(@minus,y,ybar);
    TSS=TSS*TSS';
    R2=1-diag(RSS)/diag(TSS);
    R2adj=1-(1-R2)*(T-1)/(T-k);
    Astd=reshape(sqrt(diag(results.Vml)),size(A));
    tstat=A./Astd;

    partial_r2=(tstat./sqrt(tstat.^2+T)).^2;
    
    loglik=-0.5*nvar*T*log(2*pi)-0.5*T*log(det(results.SIGml))-0.5*trace(results.SIGml\(resid*resid')); % trace(inv(results.SIGml)*resid*resid')

    crit=information_criterion();
    [normality.stat,normality.pval]=norm_test(resid');
    [autocorrelation.stat,autocorrelation.pval]=autocorr_test(resid',k,order);
    [arch.stat,arch.pval]=arch_test(resid',k,order);
    stats=struct('Astd',Astd,'tstat',tstat,'ybar',ybar,'ystd',std(y,0,2),...
        'criterion',crit,...
        'diagnostics',struct('normality',normality,...
        'autocorrelation',autocorrelation,'arch',arch),...
        'pval',2*(1-normcdf(abs(tstat))),'loglik',loglik,...
        'R2',R2,'R2adj',R2adj,'T',T,'k',k,...
        'Ftest',(T-k)/(k-1)*R2/(1-R2),'partial_r2',partial_r2);
    
end

    function crit=information_criterion()
        % maximizing the likelihood equivalent to minimizing the information
        % criterion
        crit=struct('AIC',-2*loglik/T+2*k/T,... %=log(sig2*2*pi)+1+2*k/T;
            'SC',-2*loglik/T+k*log(T)/T,...
            'HQ',-2*loglik/T+2*k*log(log(T))/T);
    end
end

function [stat,pval]=norm_test(e)
n=size(e,1);
m2=mfunc(2);
m3=mfunc(3);
m4=mfunc(4);

SK=m3/sqrt(m2.^3);
K=m4/m2.^2;
stat=n/6*SK.^2+n/24*(K-3).^2;
pval=1-chi2cdf(stat,2);
    function mj=mfunc(jj)
        mj=sum(e.^jj)/n;
    end
end

function [stat,pval]=arch_test(e,nparam,m)
[stat,pval]=autocorr_test(e.^2,nparam,m);
end

function [stat,pval]=autocorr_test(e,nparam,m)
n=numel(e);
nki=1./(n-(1:m));
r=autocorrelation(e,m);
stat=n*(n+2)*sum(nki*r.^2);
pval=1-chi2cdf(stat,n-nparam);
end

function r=autocorrelation(x,maxlag)
if nargin<2
    maxlag=1;
end
x=x(:);
r=nan(maxlag,1);
for k=1:maxlag
    r(k)=x(k+1:end)'*x(1:end-k);
end
r=r/sum(x.^2);
end
