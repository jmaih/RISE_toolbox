function [ispd,dF,iF]=CheckPositiveDefiniteness(F)

ispd=false;
iF=[];
dF=[];
if ~(any(any(isnan(F)))||any(any(isinf(F))))
    [V,D]=eig(F);
    eigvals=diag(D);
    ispd=min(eigvals)>0;
    if ispd
        dF=prod(eigvals);
        iF=V*diag(1./eigvals)*V';
    end
end

end

% all( all( M == M' ) ) & min( eig( M ) ) > 0
% Alternatively, one may use the test
%
% M = [...];     % assume M is square
% isposdef = true;
% for i=1:length(M)
%   if ( det( M(1:i, 1:i) ) <= 0 )
%     isposdef = false;
%     break;
%   end
% end
% isposdef       % 0 if false, 1 if true
