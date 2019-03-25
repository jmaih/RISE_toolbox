%  INTERNAL FUNCTION: locate_permutation locates a permutation deriving from differentiation
% 
%  ::
% 
%    ypred=locate_permutation(x,n,wrt_wise)
% 
%  Args:
% 
%     - **x** [scalar,vector|matrix]: nperm x order matrix of the permutations
%       of interest. The number of rows represents the number of permutations
%       and the number of columns the order of differentiation
%     - **n** [scalar]: number of variables in the differentiation
%     - **wrt_wise** [true|{false}]: when false, the retrieved order is that of
%       a kronecker unfolded column wise i.e. [111 112 113 121 122 123 131 132
%       133 211 212 213 221 222 223 231 232 233 311 312 313 321 322 323 331 332
%       333]. Whe true, the order is [111 211 311 121 221 321 131 231 331 112
%       212 312 122 222 322 132 232 332 113 213 313 123 223 323 133 233 333]
% 
%  Returns:
%     :
% 
%     - **ypred** [scalar|vector] : location of the permutations
% 
%