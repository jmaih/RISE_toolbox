%--- help for rsymbdiff.backdoor ---
%
%  BACKDOOR - Introduce an unknown function for differentiation
% 
%    y = backdoor(f, varargin) introduces an "alien" function f that is not
%    known to MATLAB for differentiation purposes. It creates a symbolic
%    object y that represents the alien function introduced through a backdoor.
% 
%  Input:
%    - f: The alien function to be introduced.
%    - varargin: Additional arguments that may include rsymbdiff objects.
% 
%  Output:
%    - y: A symbolic object representing the introduced alien function.
% 
%  Example:
%    f = @(x) x^2 + 1; % Define an alien function
%    y = rsymbdiff.backdoor(f, x); % Introduce the alien function to rsymbdiff
% 
%  See also: rsymbdiff, reset
%