function [pol,fval,exitflag]=orthonormal_polynomial(...
max_order,x,optimal,debug)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

a=0;
b=1;
tol=sqrt(eps);
if nargin<4
    debug=false;
    if nargin<3
        optimal=true;
        if nargin<2
            x=[];
        end
    end
end
warnstate=warning('query','all');                
warning('off')%#ok<WNOFF> %,'MATLAB:divideByZero'

options=optimset('Display','none','MaxFunEvals',100000,'tolfun',tol);
if debug
    options.Display='iter';
end
if max(max(x))>b
    error([mfilename,':: all x should be <=b'])
end
if min(min(x))<a
    error([mfilename,':: all x should be >=a'])
end
[n,N]=size(x);
nparams=sum(2:max_order+1);

if optimal && N==0
    optimal=false;
end
number_of_terms=(2:max_order+1);
%==================
if N
    POLS=ones(N,max_order,n);
    for v=1:n
        xp=x(v,:)';
        POLS(:,2,v)=xp;
        for op=2:max_order
            POLS(:,op+1,v)=POLS(:,op,v).*xp;
        end
    end
end
%==================

m=nan(n,max_order);
v=nan(n,max_order);
eta=nan(n,sum(1:max_order-1));
lb=-inf(nparams,1);
ub=inf(nparams,1);

% start at the theoretical polynomial, no matter whether there are data or
% not
nsamples=1;
pol_theory=rand(nparams,nsamples);
fval_theory=nan(nsamples,1);
for ii=1:nsamples
    [pol_theory(:,ii),fval_theory(ii),exitflag] = fmincon(@residuals,pol_theory(:,ii),[],[],[],[],lb,ub,[],options,nan,false);
end
loc=find(fval_theory==min(fval_theory),1,'first');
fval_theory=fval_theory(loc);
pol_theory=pol_theory(:,loc);

pol=pol_theory;
fval=fval_theory;
if n
    % if there are more variables, then use that polynomial as starting values
    % to find a polynomial that fits all instead of finding different
    % polynomials for different variables, which will be expensive as the
    % number of variables increases.
    if optimal
        [pol,fval_theory,exitflag] = fmincon(@all_residuals,pol_theory,[],[],[],[],lb,ub,[],options,optimal);
    end
    fval=nan(1,n);
    
    for ii=1:n
        if N>0
            [fval(ii),m(ii,:),v(ii,:),eta(ii,:)]=residuals(pol,ii,true);
        else
            if ii==1
                [~,m(ii,:),v(ii,:),eta(ii,:)]=residuals(pol,ii,false);
                fval(2:end)=fval(ones(n-1,1));
                m(2:end,:)=m(ones(n-1,1),:);
                v(2:end,:)=v(ones(n-1,1),:);
                eta(2:end,:)=eta(ones(n-1,1),:);
            end
        end
    end
end
warning(warnstate)                               

    disp('============= theory fval (target=0) =============')
    disp(fval_theory)
    disp('============= fval (target=0) =============')
    if debug
        disp(fval)
    else
        disp(['best:: ',num2str(min(fval)),' worst:: ',num2str(max(fval))])
    end
    disp('============= Means (target=0) =============')
    if debug
        disp(m)
    else
        best=abs(m)==min(min(abs(m)));
        worst=abs(m)==max(max(abs(m)));
        disp(['best:: ',num2str(min(m(best))),' worst:: ',num2str(max(m(worst)))])
    end
    disp('============= Variances (target=1) =============')
    if debug
        disp(v)
    else
        best=abs(v-1)==min(min(abs(v-1)));
        worst=abs(v-1)==max(max(abs(v-1)));
        disp(['best:: ',num2str(min(v(best))),' worst:: ',num2str(max(v(worst)))])
    end
    disp('============= Cross terms (target=0) =============')
    if debug
        disp(eta)
    else
        best=abs(eta)==min(min(abs(eta)));
        worst=abs(eta)==max(max(abs(eta)));
        disp(['best:: ',num2str(min(eta(best))),' worst:: ',num2str(max(eta(worst)))])
    end
% reshape the polynomials
% break those poly_coefs according to their order. remembering that a
% polynomial of order x has x+1 terms
pol_levels=cell(1,max_order);
offset=0;
for io=1:max_order
    pol_levels{io}=pol(offset+(1:io+1),:);
    offset=offset+io+1;
end
pol=pol_levels;
clear pol_levels

    function [res,m,v,eta]=all_residuals(params,optimal)
        res=0;
        m=nan(n,max_order);
        v=nan(n,max_order);
        eta=nan(n,sum(1:max_order-1));
        for id=1:n
            [res_i,m(id,:),v(id,:),eta(id,:)]=residuals(params,id,optimal);
            res=max(res,res_i);
        end
    end

    function [res,m,v,eta]=residuals(params,id,optimal)%
        if optimal
            [m,v,eta]=polynomial_moments(params,POLS(:,:,id));
        else
            [m,v,eta]=polynomial_moments(params,[]);
        end
        res1=sum(m.^4);
        res2=sum((v-1).^2);
        scale=1/sum(v.^4);
        res3=sum(eta.^2);
        res=scale*res1+res2+scale*res3;
        %         res=1*res1+res2+1*res3;
        res=1000*res;
    end


    function [m,v,eta]=polynomial_moments(params,xpol)
        theory=isempty(xpol);
        pol_i=nan(N,max_order);
        m=nan(1,max_order);
        v=nan(1,max_order);
        eta=nan(1,sum(1:max_order-1));
        iter=0;
        for o1=1:max_order
            offset1=sum(number_of_terms(1:o1-1));
            p1=params(offset1+(1:o1+1));
            if theory
                m(o1)=hdmr.polynomial_integration(p1,a,b);
                p11=hdmr.polynomial_multiplication(p1,p1);
                v(o1)=hdmr.polynomial_integration(p11,a,b);
                if debug==2
                    disp(['integrating for means ',num2str(m(o1)-quad(@(x)hdmr.polynomial_evaluation(p1,x),a,b,tol))])
                    disp(['integrating for variances ',num2str(v(o1)-quad(@(x)(hdmr.polynomial_evaluation(p1,x)).^2,a,b,tol))])
                end
            else
                if o1==1
                    pol_i(:,o1)=hdmr.polynomial_evaluation(p1,xpol);
                end
                m(o1)=1/N*sum(pol_i(:,o1));
                v(o1)=1/N*sum(pol_i(:,o1).^2)-m(o1)^2;
            end
            for o2=o1+1:max_order
                offset2=sum(number_of_terms(1:o2-1));
                p2=params(offset2+(1:o2+1));
                iter=iter+1;
                if theory
                    p12=hdmr.polynomial_multiplication(p1,p2);
                    eta(iter)=hdmr.polynomial_integration(p12,a,b);
                    if debug==2
                        disp(['integrating for cross terms ',num2str(eta(iter)-quad(@(x)(hdmr.polynomial_evaluation(p1,x)).*(hdmr.polynomial_evaluation(p2,x)),a,b,tol))])
                        keyboard
                    end
                else
                    if o1==1
                        pol_i(:,o2)=hdmr.polynomial_evaluation(p2,xpol);
                    end
                    eta(iter)=1/N*sum(pol_i(:,o1).*pol_i(:,o2));
                end
            end
        end
    end
end
