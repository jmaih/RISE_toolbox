function [CB] = shuffle(BC,bdim,cdim)
% perfect_shuffle produces the a perfect shuffle matrix that turns kron(B,C) into kron(C,B)
%
% Syntax
% -------
% ::
%
%   [CB] = shuffle(BC,bdim,cdim)
%
% Inputs
% -------
%
% - **BC** [matrix]: kron(B,C) of size bdim*cdim x bdim*cdim
%
% - **bdim** [numeric]: number of rows of square matrix B
%
% - **cdim** [numeric]: number of rows of square matrix C
%
% Outputs
% --------
%
% - **CB** [matrix]: kron(C,B) bdim*cdim x bdim*cdim matrix
%
% More About
% ------------
%
% In some applications, we know kron(B,C) but cannot compute kron(C,B)
% directly or it can be expensive to compute.
%
% Examples
% ---------
%
% See also:

% Reference: Carla Dee Martin (2005): "Higher-order kronecker products and
% Tensor decompositions" PhD dissertation, pp 13-14

m=bdim;
n=cdim;
mn=m*n;

[mnr,mnc]=size(BC);

if mnr~=mnc
    error('matrix must be square')
end

if mnr~=mn
    error('matrix dimensions inconsistent with bdim and cdim')
end

PAI_bc=utils.kronecker.perfect_shuffle(m,n);

CB=PAI_bc'*BC*PAI_bc;

end

