function m = nanmean(x,dim)

nans = isnan(x);
x(nans) = 0;

if nargin == 1 
    n = sum(~nans);
    n(n==0) = NaN;
    m = sum(x) ./ n;
else
    n = sum(~nans,dim);
    n(n==0) = NaN; 
    m = sum(x,dim) ./ n;
end

