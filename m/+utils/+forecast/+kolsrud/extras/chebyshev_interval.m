function ci=chebyshev_interval(y,gam,z)
if nargin<3
    z=[];
end
[N,T]=size(y);
if T>1
    c=z;
    if size(c,2)>1
        c=[];
    end
    ci=utils.forecast.kolsrud.chebyshev_band(y,gam,c);
else
    if isempty(z)
        z=utils.forecast.kolsrud.standardized_distance(y);
    end
    [~,tag]=sort(z);
    y=y(tag);
    cutoff=ceil(N*gam);
    yM=y(1:cutoff);
    ci=[
        min(yM)
        max(yM)
        ];
end

end