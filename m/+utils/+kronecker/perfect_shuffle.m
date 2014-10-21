function [Spq] = perfect_shuffle(p,q)
% perfect_shuffle produces the a perfect shuffle matrix that turns kron(B,C) into kron(C,B)
%
% Syntax
% -------
% ::
%
%   [Spq] = perfect_shuffle(p,q)
%
% Inputs
% -------
%
% - **p** [numeric]: number of rows of square matrix B
%
% - **q** [numeric]: number of rows of square matrix C
%
% Outputs
% --------
%
% - **Spq** [matrix]: p*q x p*q matrix such that
%   Spq*kron(B,C)*Spq' = kron(C,B)  
%
% More About
% ------------
%
% Examples
% ---------
%
% m1=3;n1=4; m2=5;n2=7; B=rand(m1,n1);C=rand(m2,n2);BC=kron(B,C);CB=kron(C,B);
% Sm1_m2=perfect_shuffle(m1,m2); Sn1_n2=perfect_shuffle(n1,n2); max(max(Sm1_m2*BC*Sn1_n2'-CB))==0
%
% See also:

% References: 
% - Charles F. Van Loan (2000): "The ubiquitous Kronecker product", Journal
%   of Computational and Applied Mathematics 123, pp 85-100.                 
% - Carla Dee Martin (2005): "Higher-order kronecker products and     
%   Tensor decompositions" PhD dissertation, pp 13-14

r=p*q;

Spq=nan(r);
Ir=speye(r);

offset=0;
for irow=1:q
    select_rows=irow:q:r;
    nrows=numel(select_rows);
    Spq(offset+(1:nrows),:)=Ir(select_rows,:);
    offset=offset+nrows;
end

% Spq=sparse(Spq);

end

