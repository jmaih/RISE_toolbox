function [Qx,retcode]=validate_transition_matrix(Qx)
total=sum(Qx,2);
if any(isnan(Qx(:))) || ...
        any(Qx(:)<0) || ...
        any(Qx(:)>1)|| ...
        ~(max(abs(total-1))<1e-9)
    retcode=3;
else
    retcode=0;
    Qx=bsxfun(@rdivide,Qx,total);
end
end
