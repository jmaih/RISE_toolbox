function x=TruncatedMultivariateNormalRnd(mu,SIG,lb,ub,Nsim,method,seed)
% x=TruncatedMultivariateNormalRnd(mu,SIG,lb,ub)
% x=TruncatedMultivariateNormalRnd(mu,SIG,lb,ub,Nsim)
% x=TruncatedMultivariateNormalRnd(mu,SIG,lb,ub,Nsim,'ghk')
% x=TruncatedMultivariateNormalRnd(mu,SIG,lb,ub,Nsim,'gibbs')
% x=TruncatedMultivariateNormalRnd(mu,SIG,lb,ub,Nsim,method,seed)
if nargin<7 || isempty(seed)
    seed=100*sum(clock);
    if nargin<6
        method='gibbs';
        if nargin<5
            Nsim=1000;
            if nargin<4
                error([mfilename,'::DontKnowHowToIdentifyThis'],...
                    'Number of arguments cannot be less than 5')
            end
        end
    end
elseif nargin >7
    error([mfilename,'::DontKnowHowToIdentifyThis'],...
        'Number of arguments cannot be greater than 7')
end
defaultStream = RandStream.getDefaultStream;
savedState = defaultStream.State;
s = RandStream.create('mt19937ar','seed',seed);
RandStream.setDefaultStream(s);

switch lower(method)
    case {'fast','ghk'}
        x=GewekeHajivassiliouKeane(mu,SIG,lb,ub,Nsim);
    case {'gibbs'}
        x=Gibbs(mu,SIG,lb,ub,Nsim);
    otherwise
        x=GewekeHajivassiliouKeane(mu,SIG,lb,ub,Nsim);
end

defaultStream.State = savedState;

end

function x=GewekeHajivassiliouKeane(mu,SIG,lb,ub,Nsim)
V=transpose(chol(SIG));
npar=numel(mu);
x=nan(npar,Nsim);
lb1=lb-mu;
ub1=ub-mu;
for k=1:npar
    vx=V(k,1:k-1)*x(1:k-1,:);
    tmp_low=(lb1(k)*ones(1,Nsim)-vx)/V(k,k);
    tmp_high=(ub1(k)*ones(1,Nsim)-vx)/V(k,k);
    
    x(k,:)=TruncatedRandomNormal(tmp_low,tmp_high);
end
x=mu(:,ones(1,Nsim))+V*x;

end

%====================================BEGINNING OF NEW SUB=ROUTINE ===============================================%

function x=Gibbs(mu,SIG,lb,ub,Nsim)
npar=numel(mu);
burn=100;

x0=lb+rand(npar,1).*(ub-lb);
bad=isnan(x0);
x0(bad)=mu(bad); %
V=SIG\eye(npar);

sigk=nan(npar,1);
mukfactor=nan(npar,npar-1);
for k=1:npar
    j=setdiff(1:npar,k);
    Vjj=V(j,j)-V(j,k)*V(k,j)/V(k,k); % <--- Vjj=SIG(j,j)\eye(npar-1); % huge computational savings
    mukfactor(k,:)=SIG(k,j)*Vjj;
    sigk(k)=sqrt(SIG(k,k)-SIG(k,j)*Vjj*SIG(j,k));
end

x=nan(npar,Nsim);
for isim=1:Nsim+burn
    for k=1:npar
        j=setdiff(1:npar,k);
        muk=mu(k)+mukfactor(k,:)*(x0(j)-mu(j));
        x0(k)=muk+sigk(k)*TruncatedRandomNormal((lb(k)-muk)/sigk(k),(ub(k)-muk)/sigk(k));
    end
    if isim>burn
        x(:,isim-burn)=x0;
    end
end

end

%====================================BEGINNING OF NEW SUB=ROUTINE ===============================================%

function y=TruncatedRandomNormal(lb,ub)
k_adjust=1;
crit=1e-8;

[nr,nc]=size(lb);
lb=lb(:);
ub=ub(:);
if any(isnan(lb))||any(isnan(ub))
    error([mfilename,':: nan elements'])
end
if any(lb>ub)
    error([mfilename,':: lb greater than ub'])
end
y=nan(nr*nc,1);

PHIl = NormalCumulativeDistribution(lb);
PHIr = NormalCumulativeDistribution(ub);

finite_lb_and_ub=isfinite(lb) & isfinite(ub);
same=finite_lb_and_ub &  abs(ub-lb)<crit;
tails=abs(PHIr-PHIl)<crit & ~ same; 
good_tails=tails & finite_lb_and_ub;
bad_tails=tails & ~finite_lb_and_ub;
others=~same & ~tails; clear tails

% same, no problems
if any(same)
%     disp('same')
    y(same,1)=ub(same,1);
end
% assume a uniform distribution in the tails
if any(good_tails)
%     disp('good tails')
    y(good_tails,1)=lb(good_tails,1)+(ub(good_tails,1)-lb(good_tails,1)).*rand(sum(good_tails),1);
end
% normal distribution for nice ones
if any(others)
%     disp('normal')
    y(others,1)=InverseNormalCumulativeDistribution(PHIl(others,1)+(PHIr(others,1)-PHIl(others,1)).*rand(sum(others),1));
end
% Nasty ones
if any(bad_tails)
    lbb=lb(bad_tails);
    ubb=ub(bad_tails);
    lbb(isinf(lbb))=ubb(isinf(lbb))-k_adjust;
    ubb(isinf(ubb))=lbb(isinf(ubb))+k_adjust;
    y(bad_tails,1)=lbb+(ubb-lbb).*rand(sum(bad_tails),1);
end

bad=isinf(y)|isnan(y);
if any(bad)
    disp({'lower bound','upper bound','CDF low','CDF high','draw','good tail flag','bad tail flag','other flag'
        lb(bad),ub(bad),PHIl(bad),PHIr(bad),y(bad),good_tails(bad),bad_tails(bad),others(bad)})
    error([mfilename,':: This should not happen, please contact the author of this function'])
end

y=min(max(y,lb),ub);

y=reshape(y,nr,nc);
end

%====================================BEGINNING OF NEW SUB=ROUTINE ===============================================%

function cdf=NormalCumulativeDistribution(x)
cdf=.5*(1+erf(x/sqrt(2)));
end

%====================================BEGINNING OF NEW SUB=ROUTINE ===============================================%

function x=InverseNormalCumulativeDistribution(cdf)
x=sqrt(2)*erfinv(2*cdf-1);
end

%====================================BEGINNING OF NEW SUB=ROUTINE ===============================================%
%
% function pdf=NormalDensity(x)
% pdf=1/sqrt(2*pi)*exp(-.5*x^2);
% end