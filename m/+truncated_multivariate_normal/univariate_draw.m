function y=univariate_draw(lb,ub)
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

[nr,nc]=size(lb);

lb=lb(:);

ub=ub(:);

y=nan(nr*nc,1);

PHIl = normcdf(lb);

PHIr = normcdf(ub);

same=abs(lb-ub)<10*eps(lb);

tails=~same & abs(PHIr-PHIl)<10*eps(PHIr);

good=~same & ~tails;

if any(same)
    
    y(same)=ub(same);
    
end

if any(tails)
    
    y(tails,1)=lb(tails,1)+(ub(tails,1)-ub(tails,1)).*rand(sum(tails),1);
    
end

if any(good)
    
    y(good)=InverseNormalCumulativeDistribution(PHIl(good)+(PHIr(good)-PHIl(good)).*rand(sum(good),1));
    
end

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
