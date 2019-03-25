%  INTERNAL FUNCTION: penalizes the violation of general restrictions
% 
%  ::
% 
%    p=penalize_violations(viol)
%    p=penalize_violations(viol,c)
%    p=penalize_violations(viol,c,restrictions_same_weights)
% 
%  Args:
% 
%     - **viol** [vector] : value of constraints expressed in the form g(x)<=0
%     - **c** [scalar|{10}]: penalty factor
%     - **restrictions_same_weights** [true|{false}]: decides whether all
%       constraints are equivalent or not
% 
%  Returns:
%     :
% 
%     - **p** [scalar] : penalty
% 
%  Note:
% 
%     - restrictions are expected to be of the form g(x)<=0
% 
%