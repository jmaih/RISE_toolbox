% function to minimize
fcn=@(x)sum((x(:).'-(1:numel(x))).^2);

n=5;

x0=rand(n,1);

lb=-n^2*ones(n,1);

ub=-lb;

% options
options=struct('MaxIter',1000);

[xfinal,ffinal,exitflag,H]=mc6(fcn,x0,lb,ub,options);