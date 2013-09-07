function [lpdfn,cdfn,icdfn,moments,rndfn]=dirichlet(a,b,c,d)
if nargin<4
    d=1;
    if nargin<3
        c=0;
        if nargin<2 
            if nargin==0
                if nargout>1
                    error([mfilename,':: number of output arguments cannot exceed 1 when there are no inputs'])
                end
                lpdfn=@hyperparameters;
                return
            else
                error([mfilename,':: wrong number of arguments. Must be 0,2,3 or 4'])
            end
        end
    end
end
lpdfn=@(theta)logdens_dirichlet(theta);
cdfn=@(theta)[];
icdfn=@(u)[];
moments=struct('mean',a(:)/sum(a),'sd',sqrt(a(:).*(sum(a)-a(:))/(sum(a)^2*(sum(a)+1))));
rndfn=@(n)draw_dirichlet(n);
    function lpdf=logdens_dirichlet(theta)
        % one dirichlet at a time
        if all(theta>=0)
            lpdf=gammaln(sum(a))-sum(gammaln(a(:)))+sum((a(:)-1).*log(theta(:)));
        else
            lpdf=-inf;
        end
    end
    function d=draw_dirichlet(n)
        if nargin==0
            n=1;
        end
        k=numel(a);
        u=rand(k,n);
        aa=a(:);
        d=gaminv(u,aa(:,ones(1,n)),1);
        d=bsxfun(@rdivide,d,sum(d,1));
    end
end
function [violation,space]=hyperparameters(hyper)
% -inf<a<inf, b>0
space=[1e-12,inf];
if nargin==0||isempty(hyper)
    violation=space;
else
    violation=any(hyper(:)<space(:,1))||any(hyper(:)>space(:,2));
end
end