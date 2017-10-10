function [C,R,info]=autocorr(varargin)

[C,info]=vartools.autocov(varargin{:});

stdevmat=sqrt(diag(C(:,:,1)));

bad=stdevmat<=0;

if any(bad)
    
    warning([mfilename,':: have zero or negative theoretical standard deviation'])
    
end

stdevmat=stdevmat*stdevmat';

R=C;

for ii=1:size(C,3)
    
    R(:,:,ii)=R(:,:,ii)./stdevmat;
    
end

end