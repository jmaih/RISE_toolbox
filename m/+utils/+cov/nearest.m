function vcov=nearest(vcov0,debug,farthest)

% nearest -- computes nearest covariance matrix
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

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
% quick exit
%------------
if any(oldD<too_low)
    if farthest
        oldD=abs(oldD);
    end
    D=max(oldD,too_low);
    
    vcov = V*diag(D)*V';
    
    if debug && max(abs(D-oldD))>1e-6
        warning('Covariance matrix altered')
    end
end

end