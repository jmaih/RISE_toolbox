function [A,info]=theoretical_autocorrelations(obj,ar)%,resolve_flag
if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    A=theoretical_autocovariances(obj);
    return
end
if nargin<2
    ar=[];
end

[A,info]=theoretical_autocovariances(obj,ar);%,resolve_flag

if ~info
    stdevmat=sqrt(diag(A(:,:,1)));
    bad=stdevmat<=0;
    if any(bad)
        disp({obj.varendo(bad).name}')
        error([mfilename,':: have zero or negative theoretical standard deviation'])
    end
    stdevmat=stdevmat*stdevmat';
    for ii=1:ar+1
        A(:,:,ii)=A(:,:,ii)./stdevmat;
    end
end
