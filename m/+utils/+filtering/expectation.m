%  INTERNAL FUNCTION: computes the expectation of variables from multiple regimes
% 
%  ::
% 
%    E=expectation(probs,vals)
% 
%  Args:
% 
%     - **probs** [k x T matrix]: matrix of probabilities
%     - **vals** [1 x k cell]: matrix of values for which we want to take the
%       expectation
%     - **straight** [true|{false}]: decides whether permutations should be
%       taken or not.
% 
%  Returns:
%     :
% 
%     - **E** [1 x T vector]: expecation
% 
%