function [ispd,dF,iF,F]=CheckPositiveDefiniteness(F)

if any(any(isnan(F)))
    ispd=false;
    iF=[];
    dF=[];
else
    [F,iF,dF] = hessian_conditioner(F,eps);
    ispd=isfinite(dF) && dF>0;
end


% dF=det(F);
% if dF>0
%     iF=F\eye(size(F,1));
%     if any(any(isnan(iF)))
%         ispd=false;
%     end
% else
% end

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
