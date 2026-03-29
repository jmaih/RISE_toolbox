%  Computes the initial distribution for a Markov chain.
% 
%  Args:
%    q_now_lead : Transition probabilities matrix (h x h)
%    ergodic : Flag indicating whether to use ergodic distribution (default: true)
% 
%  Returns:
%    PAI00 : Initial distribution vector (h x 1)
%    retcode : Return code (0 if successful, 308 if there are errors)
% 
%  Notes:
%    - Computes the initial distribution PAI00 for a Markov chain given the
%      transition probabilities matrix q_now_lead.
%    - If ergodic is true or if q_now_lead is diagonal, PAI00 is computed using
%      the ergodic distribution formula.
%    - Safeguards against small values in PAI00 to avoid numerical issues.
% 
%  Example:
%    q_now_lead = [0.7, 0.3; 0.4, 0.6];
%    [PAI00, retcode] = initial_markov_distribution(q_now_lead);
% 
%  See also:
%