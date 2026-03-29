%--- help for rsymbdiff.alienize ---
%
%  ALIENIZE - Introduce unknown functions through a backdoor
% 
%    eqtns = alienize(eqtns, alien_list) modifies a set of equations by
%    introducing unknown functions listed in alien_list through a backdoor
%    mechanism. The modified equations are then returned.
% 
%  Input:
%    - eqtns (char|function_handle): equations to process.
%    - alien_list (cellstr): A cell array of unknown functions to introduce.
% 
%  Output:
%    - eqtns (char): Modified equations with unknown functions introduced.
% 
%  Example:
%    eqtns = 'g(x^2 + 1)';
%    alien_list = {'g', 'h'};
%    eqtns = rsymbdiff.alienize(eqtns, alien_list);
% 
%  See also: regexprep, rsymbdiff.dealienize
%