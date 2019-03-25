%  INTERNAL FUNCTION: creates a grid of unique permutations of derivatives,
%  with non-increasing elements from left to right
% 
%  ::
% 
%    R1=derivatives_grid(R0)
% 
%  Args:
% 
%     - **R0** [vector|matrix]: initial permutation of derivatives
% 
%  Returns:
%     :
% 
%     - **R1** [matrix]: permutation indexes for derivatives of next order
% 
%  See also:
%      mygrid
% 
%