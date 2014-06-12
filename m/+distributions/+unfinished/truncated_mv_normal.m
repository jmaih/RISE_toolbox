function [theta]=truncated_mv_normal(mu,SIG,c,d,y0,n,method)
% references: 
% 1- Christian P. Robert (1995):"Simultation of truncated normal
% variables". Statistics and Computing Volume 5, Number 2 (1995), 121-125
% 2- John Geweke (): " Efficient Simulation from the Multivariate Normal
% and Student-t Distributions Subject to Linear Constraints and the
% Evaluation of Constraint Probabilities". In Computing Science and
% Statistics: Proceedings of the 23rd Symposium on the Interface, Ed. E.
% Keramidas and S. Kaufman, pp. 571-8. Fairfax Station, VA: Interface  
% Foundation of North America.

if nargin<7
    method=[];
    if nargin<6
        n=[];
        if nargin<5
            y0=[];
            if nargin<4
                d=[];
                if nargin<3
                    c=[];
                    if nargin<2
                        error([mfilename,':: at least mean and covariance should be provided'])
                    end
                end
            end
        end
    end
end

if isempty(method)
    method='robert';
end
coded_method=strcmp(method,'cdf')+2*strcmp(method,'robert')+3*strcmp(method,'geweke');
if coded_method==3 % 'geweke'
    t1=.15;
    t2=2.18;
    t3=.725;
    t4=.45;
end
mu=mu(:);
p=size(mu,1);
if isempty(n),n=1;end
if isempty(y0),y0=mu;end
if isempty(d),d=inf(p,1);end
if isempty(c),c=-inf(p,1);end

[p1,p2]=size(SIG);
if p1~=p||p2~=p
    error([mfilename,':: size of covariance matrix inconsistent with size of mean vector'])
end

theta=y0(:,ones(1,n));
V=SIG\eye(p);

sig=nan(p,1);
SIG_j_i=cell(1,p);
S_ijjj=cell(1,p);
for ii=1:p
    jj=[(1:ii-1),(ii+1:p)];
    iSIG_j_j=V(jj,jj)-V(jj,ii)*V(jj,ii)'/V(ii,ii); % <== inv(SIG(jj,jj))
    SIG_j_i{ii}=SIG(jj,ii);
    S_ijjj{ii}=SIG_j_i{ii}'*iSIG_j_j;
    sig(ii)=sqrt(SIG(ii,ii)-S_ijjj{ii}*SIG_j_i{ii});
end

iter=0;
while iter<n
    iter=iter+1;
    for ii=1:p
        jj=[(1:ii-1),(ii+1:p)];
        mu_i=mu(ii)-S_ijjj{ii}*(theta(jj,iter)-mu(jj));
        theta(ii,iter)=univariate_truncated_normal();
    end
    if iter+1<n
        theta(:,iter+1)=theta(:,iter);
    end
end

    function draw=univariate_truncated_normal()
        a1=(c(ii)-mu_i)/sig(ii);
        a2=(d(ii)-mu_i)/sig(ii);
        truncation=isfinite(a1)+2*isfinite(a2);
        switch truncation
            case 0 % unconstrained
                z=randn;
            otherwise
                if coded_method==1
                    u=rand;
                    PHI1=.5*(1+erf(a1/sqrt(2)));
                    PHI2=.5*(1+erf(a2/sqrt(2)));
                    z=sqrt(2)*erfinv(2*((PHI2-PHI1)*u+PHI1)-1);
                    if isnan(z)||isinf(z)
                        keyboard
                    end
                else
                    switch truncation
                        case 1 % lower bound is finite
                            if coded_method==2
                                z=robert_draw_left(a1);
                            else
                                if a1<=t4
                                    z=normal_rejection_sampling(a1,a2);
                                else
                                    z=exponential_rejection_sampling(a1);
                                end
                            end
                        case 2 % upper bound is finite
                            if coded_method==2
                                z=-robert_draw_left(-a2);
                            else
                                if -a2<=t4
                                    z=-normal_rejection_sampling(-a2,-a1);
                                else
                                    z=-exponential_rejection_sampling(-a2);
                                end
                            end
                        case 3 % both bounds are finite
                            if coded_method==2
                                z=robert_draw_left_right(a1,a2);
                            else
                                z=geweke_finite_bounds(a1,a2);
                            end
                    end
                end
        end
        draw=mu_i+sig(ii)*z;
    end
    function z=geweke_finite_bounds(a1,a2)
        phi_a=phi(a1);
        phi_b=phi(a2);
        if 0>=a1 && 0<=a2
            if phi_a<=t1||phi_b<=t1
                z=normal_rejection_sampling(a1,a2);
            else
                z=uniform_rejection_sampling(a1,a2);
            end
        elseif a1>0
            phi_ab=phi_a/phi_b;
            if phi_ab<=t2
                z=uniform_rejection_sampling(a1,a2);
            elseif phi_ab>t1 && a<t3
                z=half_normal_rejection_sampling(a1,a2);
            else
                z=exponential_rejection_sampling(a1);
            end
        else % if a2<0
            phi_ba=phi(-a2)/phi(-a1);
            if phi_ba<=t2
                z=-uniform_rejection_sampling(-a2,-a1);
            elseif phi_ba>t1 && -a2<t3
                z=-half_normal_rejection_sampling(-a2,-a1);
            else
                z=-exponential_rejection_sampling(-a2);
            end
        end
    end
end

function z=exponential_rejection_sampling(aa,lambda)
if nargin<2
    lambda=aa;
end
notdone=true;
while notdone
    z=exponential_draw(lambda);
    u=rand;
    if lambda<=aa
        acceptance_probability=exp(-.5*(z^2+aa^2))*exp(-lambda*z);
    else
        acceptance_probability=exp(-.5*(z-lambda)^2);
    end
    if u<=acceptance_probability
        notdone=true;
    end
end
end

function dd=exponential_draw(betta)
u=rand;
dd=-betta*log(1-u);
end

function z=uniform_rejection_sampling(a,b,phi_star)
if nargin<3
    phi_star=normpdf(0);
end
notdone=true;
while notdone
    z=a+(b-a)*rand;
    u=rand;
    if u<=phi(z)/phi_star
        notdone=false;
    end
end
end

function z=half_normal_rejection_sampling(a,b)
notdone=true;
if a<0
    error([mfilename,':: in this case, the lower bound should be positive or zero'])
end
while notdone
    z=randn;
    if z>=a && z<=b
        notdone=false;
        z=abs(z);
    end
end
end

function z=normal_rejection_sampling(a,b)
notdone=true;
while notdone
    z=randn;
    if z>=a && z<=b
        notdone=false;
    end
end
end

function z=robert_draw_left(a1)
notdone=true;
alpha_star=.5*(a1+sqrt(a1^2+4));
while notdone
    z=exponential_draw(alpha_star);
    if a1<alpha_star
        rho=exp(-.5*(alpha_star-z)^2);
    else
        rho=exp(-.5*(a1-alpha_star)^2)*exp(-.5*(alpha_star-z)^2);
    end
    if rand<=rho
        notdone=false;
    end
end
end

function z=robert_draw_left_right(a1,a2)
notdone=true;
while notdone
    z=a1+rand*(a2-a1);
    if a1<=0 && 0<=a2
        rho=exp(-.5*z^2);
    elseif a2<0
        rho=exp(-.5*(a2^2-z^2));
    elseif 0<a1
        rho=exp(-.5*(a1^2-z^2));
    end
    if rand<=rho
        notdone=false;
    end
end

end