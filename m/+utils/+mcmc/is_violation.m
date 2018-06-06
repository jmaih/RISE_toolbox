function flag=is_violation(f,penalty)
% Checks whether the fitness is too low
%
% ::
%
%   flag=is_violation(f,penalty)
%
% Args:
%
%    - **f** [numeric]: fitness
%    - **penalty** [numeric]: lowest possible fitness
%
% Returns:
%    :
%
%    - **flag** [true|false]:
%

flag=abs(f)>=penalty||~isreal(f);
if ~isreal(f)
    warning('f not real!!!')
end
end
