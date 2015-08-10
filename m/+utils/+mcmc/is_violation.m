function flag=is_violation(f,penalty)
% is_violation -- checks whether the fitness is too low
%
% Syntax
% -------
% ::
%
%   flag=is_violation(f,penalty)
%
% Inputs
% -------
%
% - **f** [numeric]: fitness
%
% - **penalty** [numeric]: lowest possible fitness 
%
% Outputs
% --------
%
% - **flag** [true|false]: 
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

flag=abs(f)>=penalty||~isreal(f);
if ~isreal(f)
    warning('f not real!!!')
end
end
