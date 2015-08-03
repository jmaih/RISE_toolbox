function [x,retcode]=sigma_points(a,P,type,scale)
if nargin<4
    scale=[];
end
if isempty(scale)
    get_scale();
end
retcode=0;
check=any(isnan(P(:)))||any(~isfinite(P(:)));% check=any(isnan(P(:)))||any(abs(P(:)>1e+5));%
if ~check
    [Pstar,check]=chol(utils.cov.project(P),'lower');
end
if check
    x=[];
    retcode=304;
else
    % symmetrizing is not enough to get a positive definite matrix
    Pstar=scale*Pstar;
    x=[bsxfun(@plus,a,Pstar),... points to the right
        bsxfun(@minus,a,Pstar)]; % points to the left
    switch type
        case 'ckf'
        case {'ukf','ddf'}
            x=[a,x]; % add center
        otherwise
            error(['unknown type of rule ',type])
    end
end

    function get_scale()
        m=size(a,1);
        switch type
            case 'ckf'
                scale=sqrt(m);
            case 'ddf'
                scale=sqrt(3);
            case 'ukf'
                k=3-m;
                scale=sqrt(m+k);
            otherwise
                error(['unknown type of rule ',type])
        end
    end
end
