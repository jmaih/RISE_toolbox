function [T,R,SS,retcode]=dsge_solve_gensys(Aplus,A0,Aminus,B,CC,T_only)
% this function solves the rational expectations model
% Aplus*X_{t+1}+A0*X_{t}+Aminus*X_{t-1}+B*E_{t}+C=0
%  The procedure adapts Chris Sims' gensys.
if nargin<6
    T_only=false;
    if nargin<5
        CC=[];
        if nargin<4
            B=[];
        end
    end
end
[g0,g1,pi,psi,c,n_expect]=form_gensys(Aplus,A0,Aminus,B,CC);

% function [G1,C,impact,fmat,fwt,ywt,gev,eu,loose]=gensys(g0,g1,c,psi,pi,div)
% System given as
%        g0*y(t)=g1*y(t-1)+c+psi*z(t)+pi*eta(t),
% with z an exogenous variable process and eta being endogenously determined
% one-step-ahead expectational errors.  Returned system is
%       y(t)=G1*y(t-1)+C+impact*z(t)+ywt*inv(I-fmat*inv(L))*fwt*z(t+1) .
% If z(t) is i.i.d., the last term drops out.
% If div is omitted from argument list, a div>1 is calculated.
% eu(1)=1 for existence, eu(2)=1 for uniqueness.  eu(1)=-1 for
% existence only with not-s.c. z; eu=[-2,-2] for coincident zeros.
% By Christopher A. Sims

T=[];SS=[];R=[];

% eu=[0;0];
realsmall=1e-6;
fixdiv=0;
n=size(g0,1);
[a b q z]=qz(g0,g1);
if ~fixdiv, div=1.01; end
nunstab=0;
zxz=0;
for i=1:n
    % ------------------div calc------------
    if ~fixdiv
        if abs(a(i,i)) > 0
            divhat=abs(b(i,i))/abs(a(i,i));
            if 1+realsmall<divhat && divhat<=div
                div=.5*(1+divhat);
            end
        end
    end
    % ----------------------------------------
    nunstab=nunstab+(abs(b(i,i))>div*abs(a(i,i)));
    if abs(a(i,i))<realsmall && abs(b(i,i))<realsmall
        zxz=1;
    end
end

if ~zxz
    [a b q z]=qzdiv(div,a,b,q,z);
end
% gev=[diag(a) diag(b)];
if zxz
    disp('Coincident zeros.  Indeterminacy and/or nonexistence.')
    %   eu=[-2;-2];
    retcode=22;
    return
end
q1=q(1:n-nunstab,:);
q2=q(n-nunstab+1:n,:);
% z1=z(:,1:n-nunstab)';
% z2=z(:,n-nunstab+1:n)';
% a2=a(n-nunstab+1:n,n-nunstab+1:n);
% b2=b(n-nunstab+1:n,n-nunstab+1:n);
etawt=q2*pi;

[ueta,deta,veta]=svd(etawt);
md=min(size(deta));
bigev=diag(deta(1:md,1:md))>realsmall; % <--- bigev=find(diag(deta(1:md,1:md))>realsmall);
ueta=ueta(:,bigev);
veta=veta(:,bigev);
deta=deta(bigev,bigev);

% eu(1) = length(bigev)>=nunstab;

etawt1 = q1 * pi;
% neta = size(pi,2);
% ndeta1 = min(n-nunstab,neta);
[ueta1,deta1,veta1]=svd(etawt1);
md=min(size(deta1));
bigev=diag(deta1(1:md,1:md))>realsmall;% <--- bigev=find(diag(deta1(1:md,1:md))>realsmall);
ueta1=ueta1(:,bigev);
veta1=veta1(:,bigev);
deta1=deta1(bigev,bigev);

if isempty(veta1)
    unique=1;
else
    retcode=21;
    loose = veta1-veta*veta'*veta1;
    [~,dl] = svd(loose);
    nloose = sum(abs(diag(dl)) > realsmall*n);
    unique = (nloose == 0);
end
if unique
    %    eu(2)=1;
    retcode=0;
end
tmat = [eye(n-nunstab) -(ueta*(deta\veta')*veta1*deta1*ueta1')'];
G0= [tmat*a; zeros(nunstab,n-nunstab) eye(nunstab)];
G1= [tmat*b; zeros(nunstab,n)];

G0I=G0\eye(n);
G1=G0I*G1;
G1=real(z*G1*z');
%====
T=G1(n_expect+1:end,n_expect+1:end);
%====
if T_only
    return
end

impact=G0I*[tmat*q*psi;zeros(nunstab,size(psi,2))];
impact=real(z*impact);
%====
R=impact(n_expect+1:end,:);
%====
if ~isempty(c)
    usix=n-nunstab+1:n;
    C=G0I*[tmat*q*c;(a(usix,usix)-b(usix,usix))\q2*c];
    C=real(z*C);
    SS=(eye(n-n_expect)-T)\C(n_expect+1:end);
end

function [G0,G1,PAI,PSI,CONS,n_expect]=form_gensys(LEAD,CURRENT,LAG,B,C)
% Recover Sims form for Gensys
% g0*y(t)=g1*y(t-1)+c+psi*z(t)+pi*eta(t)
% From Model LEAD*y_{t+1}+CURRENT*y_{t}+LAG*y_{t-1}+Bx_{t}+C=0

if nargin<3
    error([mfilename,':: Lead, current and lags should be provided at the very least'])
end
orig_n=size(CURRENT,1);
if any(size(LEAD)~=orig_n)||any(size(CURRENT)~=orig_n)||any(size(LAG)~=orig_n)
    error([mfilename,':: sizes of LEAD,CURRENT,LAG inconsistent'])
end

% find the columns with forward-looking terms
expect=any(LEAD);
n_expect=sum(expect);
n=orig_n+n_expect;
G0=zeros(n);
G0(1:orig_n,:)=[LEAD(:,expect),CURRENT];
G0(orig_n+1:end,n_expect+find(expect))=eye(n_expect);
G1=zeros(n);
G1(1:orig_n,:)=[zeros(orig_n,n_expect),-LAG];
G1(orig_n+1:end,1:n_expect)=eye(n_expect);
PAI=[zeros(orig_n,n_expect)
    eye(n_expect)];

if nargout>3
    if nargin<5
        C=[];
        if nargin<4
            B=[];
        end
    end
    if ~isempty(B)
        nshocks=size(B,2);
        PSI=[-B
            zeros(n_expect,nshocks)];
    else
        PSI=[];
    end
    if nargout>4
        if ~isempty(C)
            CONS=[C;zeros(n_expect,1)];
        else
            CONS=[];
        end
    end
end
