function [lpdfn,cdfn,icdfn,moments,rndfn]=mv_normal(a,b,c,d)
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
a=a(:);
k=numel(a);
if isequal(b,transpose(b))
    [L,p]=chol(b);
    if p==0
        C=L;
    else
        error([mfilename,':: matrix not positive definite'])
    end
    SIG=b;
else
    C=b; % assume this is the cholesky...
    SIG=b*b';
end
icdfn=@(u)[];
cdfn=@(theta)[];
lpdfn=@(theta)-k/2*log(2*pi)-1/2*log(det(SIG))-1/2*(theta(:)-a)'*(SIG\(theta(:)-a));
moments=struct('mean',a,'sd',C);
rndfn=@(n)draws_mv_normal(n);
    function d=draws_mv_normal(n)
        if nargin==0
            n=1;
        end
        d=bsxfun(@plus,a,C*randn(k,n));
    end

%=========================
    function [violation,space]=hyperparameters(~)
        % b>0
        space=[1e-12,inf];
        if nargin==0
            violation=space;
        else
            [L,p]=chol(b);
            violation=p>0;
        end
    end
end

