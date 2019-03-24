%  quantization -- discretizes a set of random variates
% 
%  ::
% 
%    [x,w] = quantization(n,xw_discrete)
% 
%  Args:
% 
%     n (vector): vector containing the number of discrete quantities in
%        each dimension of the gaussian shocks
% 
%     xw_discrete (1 x 2 cell array): The first cell is the vector of all
%        possible values of the discrete distribution. The second cell is
%        the vector of weights of each of those possible values
% 
%  Returns:
%     :
% 
%     - **x** [matrix]: grid combinations of the different variates, with
%        the number of shocks in rows and the number of combinations in
%        columns
% 
%     - **w** [vector]: probability distribution over the different
%        combinations
% 
%  Note:
% 
%     - Only works with gaussian shocks at the moment. 
% 
%