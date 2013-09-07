function Results=ordinary_least_squares(y,x,order)
if nargin<3
    order=5;
end
[T,k]=size(x);
xxi=(x'*x)\eye(k);
b=xxi*x'*y;
resid=y-x*b;
RSS=(resid'*resid);
ybar=mean(y);
TSS=(y-ybar)'*(y-ybar);
R2=1-RSS/TSS;
R2adj=1-(1-R2)*(T-1)/(T-k);
sig2=RSS/T; % maximum likelihood estimator
bstd=sqrt(sig2*diag(xxi));
tstat=b./bstd;
% partial_r2=tstat./sqrt(tstat.^2+(T-k));
partial_r2=(tstat./sqrt(tstat.^2+T)).^2;
% partial_r2_test=b;
% for ii=1:k
%     locs=[1:ii-1,ii+1:k];
%     xi=x(:,locs);
%     res_y=y-xi*(xi\y);
%     res_x=x(:,ii)-xi*(xi\x(:,ii));
%     partial_r2_test(ii)=corr(res_y,res_x)^2;
% end
% disp(max(abs(partial_r2-partial_r2_test)))
% keyboard
loglik=-.5*T*(log(sig2*2*pi)+1); % <--- loglik=-.5*T*log(sig2*2*pi)-0.5/sig2*sum(resid.^2);
crit=information_criterion();
[normality.stat,normality.pval]=norm_test(resid);
[autocorrelation.stat,autocorrelation.pval]=autocorr_test(resid,k,order);
[arch.stat,arch.pval]=arch_test(resid,k,order);
Results=struct('b',b,'e',resid,'sig2',sig2,'bstd',bstd,'tstat',tstat,...
    'criterion',crit,...
    'diagnostics',struct('normality',normality,...
    'autocorrelation',autocorrelation,'arch',arch),...
    'pval',2*(1-normcdf(abs(tstat))),'loglik',loglik,'RSS',RSS,...
    'ybar',ybar,'ystd',std(y),'R2',R2,'R2adj',R2adj,'T',T,'k',k,...
    'Ftest',(T-k)/(k-1)*R2/(1-R2),'partial_r2',partial_r2);

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
