function p=mypermutation(v,p)
if nargin<2
    p={[]};
end
n=numel(v);
np=numel(p);
if n==1
    for ip=1:np
        p{ip}=[p{ip},v];
    end
else
    iter=0;
    for ii=1:n
        vi=v;
        vi(ii)=[];
        ppi=mypermutation(vi,p);
        if ii==1
            nppi=numel(ppi);
            pp_nbr=nppi*n;
            pp=cell(pp_nbr,1);
        end
        for jj=1:nppi
            iter=iter+1;
            pp{iter}=[v(ii),ppi{jj}];
        end
    end
    p=pp;
end

end