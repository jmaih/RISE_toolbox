function flag=is_computable(x)
% is_computable -- checks whether the input is valid for computation
%
% Syntax
% -------
% ::
%
%   flag=is_computable(x)
%
% Inputs
% -------
%
% - **x** [empty|matrix]: matrix to be checked for non-emptiness and for
% any non-zeros
%
% Outputs
% --------
%
% - **flag** [true|false]: true if the matrix is not empty and has at least
% one non-zero element
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:



flag=~isempty(x) && any(x(:)~=0);

end