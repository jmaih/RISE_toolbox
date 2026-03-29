%--- help for rsymbdiff.dealienize ---
%
%  DEALIENIZE - Remove the backdoor mechanism around unknown functions
% 
%    eqtns = dealienize(eqtns, alien_list) 
% 
%  Input:
%    - eqtns (char|function_handle): equations to process.
%    - alien_list (cellstr): A cell array of unknown functions to introduce.
% 
%  Output:
%    - eqtns (char): Modified equations with unknown functions introduced.
% 
%  See also: regexprep, rsymbdiff.alienize
%