%--- help for build_model ---
%
%  Create inputs for model setting and estimation of a structural VAR identified by Choleski restrictions.
% 
%  Args:
%     - m (string): nvmc or nvABCmc, where n is the number of states for the
%       variances, m the number of states for the coefficients. A, B, C are the
%       equations identified by variable names.
%     - vnames (cell of string): list of the endogenous variables ordered according
%       to their position in the VAR
%     - options (struct): [{}] default options for beta and dirichlet
%       distributions
% 
%  Returns:
%     :
%        - mc (struct): structure containing information about the markov chains
%        - switch_prior (struct): structure containing the priors on the transition
%          probabilities
%        - restr (cell array): list of (Choleski) restrictions. state-identifying
%          restrictions are not added yet
% 
%