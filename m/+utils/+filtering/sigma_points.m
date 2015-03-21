function [x]=sigma_points(a,P,type,scale)
if nargin<4
    scale=[];
end
if isempty(scale)
    get_scale();
end
% symmetrizing is not enough to get a positive definite matrix
Pstar=scale*chol(utils.cov.nearest(P),'lower');
x=[bsxfun(@plus,a,Pstar),... points to the right
    bsxfun(@minus,a,Pstar)]; % points to the left
switch type
    case 'ckf'
    case {'ukf','ddf'}
        x=[a,x]; % add center
    otherwise
        error(['unknown type of rule ',type])
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
