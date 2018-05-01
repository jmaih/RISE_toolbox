function [f,xx]=empirical_cdf(x,lb,ub,N)
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