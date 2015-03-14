function vcov=nearest(vcov0,debug,farthest)

% nearest -- computes nearest covariance matrix
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

if nargin<3
    farthest=false;
    if nargin<2
        debug=false;
    end
end

vcov=.5*(vcov0+vcov0.');

[V,D] = eig(vcov);

oldD=diag(D);

% quick exit
%------------
if any(oldD<=0)
    if farthest
        oldD=abs(oldD);
    end
    D=max(oldD,sqrt(eps));
    
    vcov = V*diag(max(D,sqrt(eps)))*V';
    
    if debug && max(abs(D-oldD))>1e-6
        warning('Covariance matrix altered')
    end
end

end