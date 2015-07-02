function C=kronall(varargin)
% kronall -- multikronecker
%
% Syntax
% -------
% ::
%
%   C=kronall(A1,A2,...,An)
%
% Inputs
% -------
%
% - **Ai** [matrix]: input for the kronecker product
%
% Outputs
% --------
%
% - **C** [matrix]: kron(A1,kron(A2,kron(A3,...)))
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: tensorperm

test=true;
nargs=length(varargin);
if test
    C=varargin{end};
    for ii=nargs-1:-1:1
        C=kron(varargin{ii},C);
    end
else
    C=varargin{1};
    for ii=2:nargs
        C=kron(C,varargin{ii});
    end
end


end