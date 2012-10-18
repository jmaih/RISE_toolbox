function H = hessian_conditioner(H0,tol,loop)
% based on ddom.m
% returns a diagonally dominant matrix H by modifying the
% diagonal of H0.

if nargin<3
    loop=false;
    if nargin<2
        tol=100*eps;
    end
end
[m,n] = size(H0);
if m~=n
    error([mfilename,':: matrix must be square'])
end

H=H0;
if loop
    for ii=1:n
        d=H0(ii,ii);
        a=abs(d);
        f=0;
        for jj=1:n
            if ii~=jj
                f=f+abs(H0(ii,jj));
            end
        end
        if f>=a
            aii=(1+tol)*max(f,tol);
            if d<0
                aii=-aii;
            end
            H(ii,ii)=aii;
        end
    end
else
    d = diag(H0);
    a = abs(d);
    f = sum(abs(H0),2)-a;
    ii = find(f >= a);
    k = ii + (ii-1)*m;
    s = 2*(d(ii)>=0)-1;
    H(k) = (1+tol)*s.*max(f(ii),tol);
end
