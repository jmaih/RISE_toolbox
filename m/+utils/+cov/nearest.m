function vcov=nearest(vcov0)
% nearest covariance matrix
vcov=.5*(vcov0+vcov0.');
[V,D] = eig(vcov);
vcov = V*diag(max(diag(D),sqrt(eps)))*V';
if max(abs(vcov(:)-vcov0(:)))>1e-6
    warning('Covariance matrix altered')
end
end