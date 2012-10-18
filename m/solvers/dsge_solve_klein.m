function [T,R,SS,retcode]=dsge_solve_klein(Aplus,A0,Aminus,B,junk,T_only)
% this function solves the rational expectations model
% Aplus*X_{t+1}+A0*X_{t}+Aminus*X_{t-1}+B*E_{t}+C=0 using Klein's method.
% The matrices are stacked to yield the Klein/Gensys form
if nargin<6
    T_only=false;
    if nargin<5
        if nargin<4
            B=[];
        end
    end
end
% set some options
klein_order=false; % no need to order, really
direct=false; % uses bigger matrices...

% initialize output
R=[];
SS=[];

n=size(A0,2);

if klein_order
    static=~any(Aplus) & ~any(Aminus);
    pred=~any(Aplus) & any(Aminus);
    both=any(Aminus) & any(Aplus);
    fl=any(Aplus) & ~any(Aminus);
    order=[find(static),find(pred),find(both),find(fl)];
else
    order=1:n;
end
[AAplus,AA0,AAminus]=reorder(order,Aplus,A0,Aminus);
inv_order(order)=1:n;

if direct
    [D,E]=form_klein_full(AAplus,AA0,AAminus,B);
    nshocks=size(B,2);
    nstates=nshocks+size(AA0,2);
else
    [D,E,F]=form_klein(AAplus,AA0,AAminus,B);
    nstates=size(AA0,2);
end

[T,junk,retcode]=solve_klein(D,E,nstates);

if ~retcode
    if direct
        T0=T;
        T=T0(nshocks+1:end,nshocks+1:end);
    end
   % reorder rows and columns
    T=T(inv_order,inv_order);
    if ~T_only
        if direct
            R=T0(nshocks+1:end,1:nshocks);
        else
            R=(AAplus*T+AA0)\B;
        end
        % reorder rows
        R=R(inv_order,:);
    end
end


function varargout=reorder(order,varargin)
for ii=1:numel(varargin)
    varargout{ii}=varargin{ii}(:,order);
end

function [sol_T,eigval,retcode]=solve_klein(AA,BB,nb,tol)
if nargin<4
    tol=sqrt(eps);
end
% AA*X_{t+1}+BB*X_{t+1}=0;
[SS,TT,QQ,ZZ] = qz(AA,BB,'real');
% Ordered inverse eigvals.
eigval = ordeig(SS,TT); % eigval = -ordeig(SS,TT);
% eigval = eigval(:).';
stable = abs(eigval) >= 1 + tol;
unit = abs(abs(eigval)-1) < tol;
% Clusters of unit, stable, and unstable eigenvalues.
clusters = zeros(size(eigval));
% Unit roots first.
clusters(unit) = 2;
% Stable roots second.
clusters(stable) = 1;
% Unstable roots last.
% Re-order by the clusters.
[SS,TT,junk,ZZ] = ordqz(SS,TT,QQ,ZZ,clusters);

% Undo the eigval inversion.
infeigval = eigval == 0;
eigval(~infeigval) = 1./eigval(~infeigval);
eigval(infeigval) = Inf;
nunit = sum(unit);
nstable = sum(stable);
%%

% Check BK saddle-path condition.
% f=[];
sol_T=[];
if nb == nstable + nunit
    retcode=0; %disp('unique solution')
    S11 = SS(1:nb,1:nb);
    %     S12 = SS(1:nb,nb+1:end);
    %     S22 = SS(nb+1:end,nb+1:end);
    T11 = TT(1:nb,1:nb);
    %     T12 = TT(1:nb,nb+1:end);
    %     T22 = TT(nb+1:end,nb+1:end);
    Z11 = ZZ(1:nb,1:nb);
    %     Z12 = ZZ(1:nb,nb+1:end);
    %     Z21 = ZZ(nb+1:end,1:nb);
    %     Z22 = ZZ(nb+1:end,nb+1:end);
    Z11i=Z11\eye(nb); 
    dyn = -(S11\T11); % mind the minus sign
    %     f = Z21*Z11i; % explosive part, not needed
    sol_T = Z11*dyn*Z11i;
elseif nb > nstable + nunit
    retcode=22;
    %     disp('no solution')
else
    %     disp('multiple solutions')
    retcode=21;
end


function [D,E,F,n_expect]=form_klein(LEAD,CURRENT,LAG,B)
% Recover Klein form for Solab
% D*y(t)+E*y(t-1)+F*eta(t)+G=0
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
nshocks=size(B,2);
if nargout
    D=zeros(n);
    D(1:orig_n,:)=[CURRENT,LEAD(:,expect)];
    D(orig_n+1:end,find(expect))=eye(n_expect); %#ok<FNDSB>
    if nargout>1
        E=zeros(n);
        E(1:orig_n,1:orig_n)=LAG;
        E(orig_n+1:end,orig_n+1:end)=-eye(n_expect);
        if nargout>2
            F=B;
            if ~isempty(B)
                F=[F
                    zeros(n_expect,nshocks)];
            end
            if nargout>3
                error([mfilename,':: number of arguments cannot exceed 3'])
            end
        end
    end
end

function [D,E,n_expect]=form_klein_full(LEAD,CURRENT,LAG,B)
% Recover Klein form for Solab
% D*y(t)+E*y(t-1)+F*eta(t)+G=0
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
nshocks=size(B,2);
n=nshocks+orig_n+n_expect;
if nargout
    D=zeros(n);
    D(1:orig_n,nshocks+1:end)=[CURRENT,LEAD(:,expect)];
    D(orig_n+(1:n_expect),nshocks+find(expect))=eye(n_expect); 
    D(orig_n+n_expect+1:end,1:nshocks)=eye(nshocks); 
    if nargout>1
        E=zeros(n);
        E(1:orig_n,1:nshocks+orig_n)=[B,LAG];
        E(orig_n+(1:n_expect),nshocks+orig_n+1:end)=-eye(n_expect);
        if nargout>3
            error([mfilename,':: number of arguments cannot exceed 3'])
        end
    end
end

% % % % % reduce=false;
% % % % % if reduce
% % % % %     [qq,rr]=qr(AA0);
% % % % %     nstat=sum(static);
% % % % %     npred=sum(pred);
% % % % %     nboth=sum(both);
% % % % %     nfl=sum(fl);
% % % % %     qAplus=qq'*AAplus;
% % % % %     qA0=qq'*AA0;
% % % % %     qAminus=qq'*AAminus;
% % % % %     qB=qq'*B;
% % % % %     if isempty(C)
% % % % %         C=zeros(n,1);
% % % % %     end
% % % % %     qC=qq'*C;
% % % % %     
% % % % %     rAplus=qAplus(nstat+1:end,nstat+1:end);
% % % % %     rA0=qA0(nstat+1:end,nstat+1:end);
% % % % %     rAminus=qAminus(nstat+1:end,nstat+1:end);
% % % % %     rB=qB(nstat+1:end,:);
% % % % %     rC=qC(nstat+1:end,:);
% % % % %     [D,E,F,CONS]=form_klein(rAplus,rA0,rAminus,rB,rC);
% % % % %     nstates=size(rA0,1);
% % % % % end

% % % % function [f,p,lambda,retcode] = solab(a,b,nstates)
% % % % % Adapted from original solab by Paul Klein. Replaces the call to reorder
% % % % % by a call to ordqz
% % % % % Purpose: Solves for the recursive representation of the stable solution
% % % % % to a system of linear difference equations.
% % % % % Inputs: Two square matrices a and b and a natural number nstates
% % % % % a and b are the coefficient matrices of the difference equation
% % % % % a*x(t+1)+b*x(t)= 0
% % % % % where x(t) is arranged so that the state variables come first,
% % % % % nstates is the number of state variables.
% % % % % Outputs: the decision rule f and the law of motion p. If we write
% % % % % x(t) = [k(t);u(t)] where k(t) contains precisely the state bles, then
% % % % % u(t) = f*k(t) and
% % % % % k(t+1) = p*k(t).
% % % % [s,t,q,z] = qz(a,b,'real'); % upper triangular factorization
% % % % % % pred=any(b);
% % % % % % nstates=sum(pred);
% % % % % of the matrix pencil b-za
% % % % % [s,t,~,z] = reorder(s,t,q,z); % reordering of generalized eigenvalues in
% % % % n=size(a,1);
% % % % dt=diag(t);
% % % % ds=diag(s);
% % % % lambda_inv=ds./dt;
% % % % lambda_inv(dt==0)=inf*sign(ds(dt==0));
% % % % tol=1e-7;
% % % % clusters=zeros(1,n);
% % % % % unit roots
% % % % clusters(abs(abs(lambda_inv)-1)<tol)=2;
% % % % % stationary
% % % % clusters(abs(lambda_inv)>=1+tol)=1;
% % % % % the rest... explosive
% % % % % clusters(abs(lambda)>1+tol)=0;
% % % % [ss,tt,~,zz] = ordqz(s,t,q,z,clusters);
% % % % % reordering of generalized eigenvalues like Benes, unit roots,stable and
% % % % % then explosive
% % % % lambda = ordeig(ss,tt);
% % % % % % check that things are ordered as expected
% % % % % [lambda,abs(ordeig(ss,tt))]
% % % % z21 = zz(nstates+1:end,1:nstates);
% % % % z11 = zz(1:nstates,1:nstates);
% % % % 
% % % % nexp=size(a,1)-nstates;
% % % % 
% % % % retcode=0;
% % % % f=[];
% % % % p=[];
% % % % if rank(z11)<nstates % disp('Invertibility condition violated')
% % % %     retcode =1;
% % % % elseif abs(tt(nstates,nstates))>abs(s(nstates,nstates)) ||...
% % % %         (nexp && abs(tt(nstates+1,nstates+1))>abs(s(nstates+1,nstates+1)))
% % % %     retcode =2; %'Wrong number of stable eigenvalues.'
% % % % else
% % % %     z11i = z11\eye(nstates);
% % % %     s11 = ss(1:nstates,1:nstates);
% % % %     t11 = tt(1:nstates,1:nstates);
% % % %     dyn = s11\t11;
% % % %     f = real(z21*z11i);
% % % %     p = -real(z11*dyn*z11i);
% % % % end
