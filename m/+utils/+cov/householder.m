function S = householder(A)

% Householder - finds S from A such that S*S' = A*A' with S square and A not
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
%     S=triag(A) uses a Householder transformation of the rectangular
%     matrix A to produce a square and upper triangular matrix S with
%     the property S*S' = A*A'.
%
% Examples
% ---------
%
% See also: 

% HOUSEHOLDER
% -----

%     Literature:
%     The algorithm is described in Grewal & Andrews: "Kalman Filtering - 
%     Theory and Practice" Prentice-Hall, 1993, pp. 234. However, there 
%     are certain cases in which their implementation fails. These have been
%     handled as described in G.W. Stewart: "Matrix Algorithms, 
%     Vol. 1: Basic Decompositions", SIAM 1988.
%
% Original code has name triag.m as was written by: Magnus Norgaard,
% IMM/IAU, Technical University of Denmark 

[n,rn]=size(A);          % Rows and columns of A
r = rn-n;
v = zeros(n,1);
u = zeros(1,rn);
S = zeros(n,n);

for k=n:-1:1,
   u(1:r+k) = A(k,1:r+k);
   ny = norm(u(1:r+k));
   if ny==0,
     u(r+k) = sqrt(2);
   else
      u(1:r+k) = u(1:r+k)/ny;
      if u(r+k)>=0,  
        u(r+k) = u(r+k) + 1;
        ny = -ny;
      else
        u(r+k) = u(r+k) - 1;
      end
      u(1:r+k) = u(1:r+k)/sqrt(abs(u(r+k)));
   end
   v(1:k-1,1) = A(1:k-1,1:r+k)*u(1,1:r+k)';
   A(1:k-1,1:r+k) = A(1:k-1,1:r+k) - v(1:k-1,1)*u(1,1:r+k);
   S(1:k-1,k) = A(1:k-1,r+k);
   S(k,k)=ny;
end