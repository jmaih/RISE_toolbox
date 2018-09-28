function [vcov,is_altered]=nearest(vcov0,debug,farthest)
% INTERNAL FUNCTION: Computes nearest covariance matrix
%

if nargin<3
    
    farthest=false;
    
    if nargin<2
        
        debug=false;
        
    end
    
end

vcov=.5*(vcov0+vcov0.');

[V,D] = eig(vcov);

oldD=diag(D);

too_low=sqrt(eps);

is_altered=false;
% quick exit
%------------
if any(oldD<too_low)
    
    if farthest
        
        oldD=abs(oldD);
        
    end
    
    D=max(oldD,too_low);

    vcov = V*diag(D)*V';
    
    is_altered=max(abs(D-oldD))>1e-6;

    if debug && is_altered
        
        warning('Covariance matrix altered')
        
    end
    
end

end