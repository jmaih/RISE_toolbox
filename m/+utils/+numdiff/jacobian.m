function J=jacobian(func,x,varargin)
tol=eps^(1/3);

x=x(:);
h=tol*max(1,x);
xp=x+h;
xm=x-h;
h=xp-xm;
n=numel(x);

x=x(:,ones(n,1));

diag_terms=(0:n-1)*n+(1:n);
xxp=x;
xxp(diag_terms)=xp;
xxm=x;
xxm(diag_terms)=xm;

J=func(xxp,varargin{:})-func(xxm,varargin{:});
nrows=size(J,1);
hp=h';
J=J./hp(ones(nrows,1),:);

end

