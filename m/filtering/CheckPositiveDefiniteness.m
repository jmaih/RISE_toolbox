function [ispd,dF,iF]=CheckPositiveDefiniteness(F)
ispd=true;
dF=det(F);
if dF>0
    iF=F\eye(size(F,1));
    if any(any(isnan(iF)))
        ispd=false;
    end
else
	iF=[];
    dF=[];
    ispd=false;
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
