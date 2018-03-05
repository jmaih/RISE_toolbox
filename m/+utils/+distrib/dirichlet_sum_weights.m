function sum_aij=dirichlet_sum_weights(h)
% DIRICHLET_SUM_WEIGHTS -- computes the sum of the "weights" of the
% transformed/inverse dirichlet
%
% Syntax
% -------
% ::
%
%   sum_aij=DIRICHLET_SUM_WEIGHTS(h)
%
% Inputs
% -------
%
% - **h** [integer]: number of regimes
%
% Outputs
% --------
%
% - **sum_aij** [numeric]: sum of weights
%
% More About
% ------------
%
% - Under estimation then, the transformed parameters will be
% xj=aj/sum_aij, where sum_aij=a_ii+sum(aj). Hence x_ii=1-sum(xj). It is
% well possible that sum(xj)=1, in which case x_ii=0
%
% Examples
% ---------
%
% See also: 

sum_aij=1*(h-1);

end