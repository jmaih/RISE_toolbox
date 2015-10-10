function cb=chebyshev_band(y,gam,c)
if nargin<3
    c=utils.forecast.kolsrud.chebyshev_distance(y);
end
N=numel(c);
[~,tag]=sort(c);
y=y(tag,:);
cutoff=ceil(N*gam);
yM=y(1:cutoff,:);
cb=[
    min(yM,[],1)
    max(yM,[],1)
    ];
end