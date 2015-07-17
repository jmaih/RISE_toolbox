function X = posdef(A,tol)


[U,S,V] = svd(A,'econ');
s = diag(S);
if nargin < 2 
    tol = max(size(A)) * eps(norm(s,inf));
end
r1 = sum(s > tol)+1;
V(:,r1:end) = [];
U(:,r1:end) = [];
s(r1:end) = [];
X = bsxfun(@times,V,s(:).')*U';
