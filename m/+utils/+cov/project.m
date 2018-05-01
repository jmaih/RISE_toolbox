function vcov = project(vcov0,e_min,e_max)
% project -- projection of covariance matrix such the eigenvalues are
% sandwiched.
%
% ::
%
%
% Args:
%
%    - **vcov0** [matrix]: initial covariance matrix
%
%    - **e_min** [[]|{sqrt(eps)}]: scalar such that the minimum eigenvalue of vcov
%    is greater than or equal to "e_min".
%
%    - **e_max** [[]|{1/e_min}]: scalar such that maximum eigenvalue of vcov
%    is less than or equal to "e_max"
%
% Returns:
%    :
%
%    - **vcov** [matrix]: updated covariance matrix
%
% Note:
%
% Example:
%
%    See also:

if nargin<3
    e_max=[];
    if nargin<2
        e_min=[];
    end
end

if isempty(e_min)
    e_min=sqrt(eps);
end
if isempty(e_max)
    e_max=1/e_min;
end

vcov=.5*(vcov0+vcov0.');

[V,D] = eig(vcov);

oldD=diag(D);

% quick exit
%------------
if any(oldD<e_min)||any(oldD>e_max)
    D=max(oldD,e_min);
    D=min(D,e_max);
    vcov = V*diag(D)*V';
end


end
