function [f,xx]=empirical_cdf(x,lb,ub,N)
% INTERNAL FUNCTION: Approximate empirical cumulative distribution function
%
% ::
%
%    [f, xx] = empirical_cdf(x,lu,ub,N);
%    [f, xx] = empirical_cdf(x,lu,ub);
%    [f, xx] = empirical_cdf(x,lu);
%    [f, xx] = empirical_cdf(x);
%
% Args:
%    x (vector of double): data points to create empirical CDF of
%    lb (double): lower bound of the range (default min(x))
%    ub (double): upper bound of the range (default max(x))
%    N (integer): number of bins/grid points to use to approximate the empirical CDF (default 250)
%
% Returns:
%    :
%    xx (N x 1 double): knot points of the emprical CDF (uniformly distributed between lb and ub)
%    f (N x 1 double): value of the CDF at the given point
%

if nargin<4
    N=[];
	if nargin<3
		ub=[];
		if nargin<2
			lb=[];
		end
	end
end
if isempty(N)
	N=250;
end
if isempty(ub)
	ub=max(x);
end
if isempty(lb)
	lb=min(x);
end
npar=numel(x);
xx=transpose(linspace(lb,ub,N));
f=zeros(N,1);
for i=1:N
    if ~isempty(x)
        target=x<=xx(i);
        x=x(~target);
        f(i)=sum(target);
    end
    if i>1
        f(i)=f(i)+f(i-1);
    end
end

f=f/npar;
