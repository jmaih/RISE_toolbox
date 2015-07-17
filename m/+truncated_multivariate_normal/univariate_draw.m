function y=univariate_draw(lb,ub)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

[nr,nc]=size(lb);
lb=lb(:);
ub=ub(:);
y=nan(nr*nc,1);
PHIl = normcdf(lb);
PHIr = normcdf(ub);
same=abs(lb-ub)<10*eps(lb);
tails=~same & abs(PHIr-PHIl)<10*eps(PHIr);
good=~same & ~tails;

y(same)=ub(same);
y(tails)=lb(tails)+(ub(tails)-ub(tails)).*rand(sum(tails),1);
y(good)=InverseNormalCumulativeDistribution(PHIl(good)+(PHIr(good)-PHIl(good)).*rand(sum(good),1));

y=reshape(y,nr,nc);
end

%====================================BEGINNING OF NEW SUB=ROUTINE ===============================================%

function cdf=normcdf(x)
cdf=.5*(1+erf(x/sqrt(2)));
end

%====================================BEGINNING OF NEW SUB=ROUTINE ===============================================%

function x=InverseNormalCumulativeDistribution(cdf)
x=sqrt(2)*erfinv(2*cdf-1);
end
