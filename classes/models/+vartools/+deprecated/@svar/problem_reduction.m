function [T,Q,n,h]=problem_reduction(obj)
% PROBLEM_REDUCTION - prepares the solution for checking Mean Square
% Stability of the system
%
% Syntax
% -------
% ::
%
%   [T,Q,n,h]=problem_reduction(obj)
%
% Inputs
% -------
%
% - **obj** [rfvar|svar]: solved model object
%
% Outputs
% --------
%
% - **T** [1 x h cell]: Companion form for different regimes. Each cell
% contains an n x n matrix
%
% - **Q** [h x h matrix]: Transition matrix
%
% - **n** [integer]: number of endogenous variable
%
% - **h** [integer]: number of regimes
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: DSGE/PROBLEM_REDUCTION

T=load_solution(obj,'ov',false);

Q=obj.solution.transition_matrices.Q;

h=numel(T);

n=size(T{1},1);

for ireg=1:h
    T{ireg}=T{ireg}(:,1:n);
end

end