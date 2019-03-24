%  Class for Monte Carlo Filtering
% 
%  Attributes:
% 
%     - lb : lower bounds for the parameters
%     - ub : upper bounds for the parameters
%     - nsim : number of simulations
%     - procedure : sampling procedure [ {'uniform'} | 'latin_hypercube' | 'sobol' | 'halton' | user-defined]
%     - parameter_names : names of the parameters
%     - samples : parameter draws
%     - is_behaved : boolean flag for behaved parameter vectors
%     - nparam : number of parameters
%     - is_sampled : true if draws are available
%     - check_behavior : checks whether the vectors should be true or false
%     - number_of_outputs : number of outputs to check_behavior
%     - user_outputs : sampled additional outputs
%     - known_procedures : known sampling procedures
% 
% 
%
%    Reference page in Doc Center
%       doc mcf
%
%