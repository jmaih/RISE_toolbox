function [PAI_bc] = perfect_shuffle(bdim,cdim)
% perfect_shuffle produces the a perfect shuffle matrix that turns kron(B,C) into kron(C,B)
%
% Syntax
% -------
% ::
%
%   [PAI_bc] = perfect_shuffle(bdim,cdim)
%
% Inputs
% -------
%
% - **bdim** [numeric]: number of rows of square matrix B
%
% - **cdim** [numeric]: number of rows of square matrix C
%
% Outputs
% --------
%
% - **PAI_bc** [matrix]: bdim*cdim x bdim*cdim matrix such that
%   PAI_bc'*kron(B,C)*PAI_bc = kron(C,B)  
%
% More About
% ------------
%
% Examples
% ---------
%
% nb=5;nc=7; B=rand(nb);C=rand(nc);BC=kron(B,C);CB=kron(C,B);
% PAI_bc=perfect_shuffle(nb,nc); max(max(PAI_bc'*BC*PAI_bc-CB))==0
%
% See also:

% Reference: Carla Dee Martin (2005): "Higher-order kronecker products and
% Tensor decompositions" PhD dissertation, pp 13-14

m=bdim;
n=cdim;

mn=m*n;

PAI_bc=nan(mn);
Imn=eye(mn);

offset=0;
for irow=1:m
    select_rows=irow:m:mn;
    nrows=numel(select_rows);
    PAI_bc(offset+(1:nrows),:)=Imn(select_rows,:);
    offset=offset+nrows;
end

PAI_bc=sparse(PAI_bc);

end

