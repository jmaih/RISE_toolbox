%  INTERNAL FUNCTION: computes the sum of the "weights" of the transformed/inverse dirichlet
% 
%  ::
% 
%    sum_aij=DIRICHLET_SUM_WEIGHTS(h)
% 
%  Args:
% 
%     - **h** [integer]: number of regimes
% 
%  Returns:
%     :
% 
%     - **sum_aij** [numeric]: sum of weights
% 
%  Note:
% 
%     - Under estimation then, the transformed parameters will be
%       xj=aj/sum_aij, where sum_aij=a_ii+sum(aj). Hence x_ii=1-sum(xj). It is
%       well possible that sum(xj)=1, in which case x_ii=0
% 
%