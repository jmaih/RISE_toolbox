%  INTERNAL FUNCTION: computes the regimes from one or multiple markov chains
% 
%  ::
% 
%    Regimes=chain_grid(a); : one chain with "a" states
% 
%    Regimes=chain_grid([a,b]); : two chains with "a" states for the first
%    and "b" states for the second
% 
%    Regimes=chain_grid([a1,a2,...,an]); : n chains with "a1" states for the
%    first, "a2" states for the second,...,"an" states for the nth
% 
%    [Regimes,Journal]=chain_grid(...);
% 
%  Args:
% 
%     - **a** [scalar|vector]: number of states in each one of the independent
%       Markov chains
% 
%  Returns:
%     :
% 
%     - **Regimes** [matrix]: Combinations of all states. The matrix is of size
%       prod([a1,a2,...,an]) x n
%     - **Journal** [cell array]: Description of the transition matrix. Each
%       cell contains a matrix in which the first row is the regime today and
%       the second one is the regime tomorrow. The cell array is of size
%       prod([a1,a2,...,an]) x prod([a1,a2,...,an])
% 
%  Example:
% 
%     ::
% 
%        Regimes=chain_grid(10); % one chain with 10 regimes
%        Regimes=chain_grid([10,5]); % 2 chains with 5 regimes for the second one
%        Regimes=chain_grid([10,5,3]); % 3 chains with 3 regimes on the third one
% 
%