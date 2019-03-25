%  INTERNAL FUNCTION: Finds S from A such that S*S' = A*A' with S square and A not
% 
%  Note:
%     S=triag(A) uses a Householder transformation of the rectangular
%     matrix A to produce a square and upper triangular matrix S with
%     the property S*S' = A*A'.
% 
%  References:
%      - The algorithm is described in Grewal & Andrews: "Kalman Filtering -
%        Theory and Practice" Prentice-Hall, 1993, pp. 234. However, there
%        are certain cases in which their implementation fails. These have been
%        handled as described in G.W. Stewart: "Matrix Algorithms,
%        Vol. 1: Basic Decompositions", SIAM 1988.
%      - Original code has name triag.m as was written by: Magnus Norgaard,
%        IMM/IAU, Technical University of Denmark
% 
%