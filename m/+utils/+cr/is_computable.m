function flag=is_computable(varargin)
% is_computable -- checks whether the input is valid for computation
%
% Syntax
% -------
% ::
%
%   flag=is_computable(x1)
%
%   flag=is_computable(x1,x2,...xn)
%
% Inputs
% -------
%
% - **x** [empty|matrix]: matrix to be checked for non-emptiness and for
% any non-zeros. If there are many x's, 
%
% Outputs
% --------
%
% - **flag** [true|false]: true if all the matrices are not empty and have
% at least one non-zero element
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

flag=true;
ivar=0;
nargs=length(varargin);
while flag && ivar < nargs
    ivar=ivar+1;
    x=varargin{ivar};
    flag=~isempty(x) && any(x(:)~=0);
end

end