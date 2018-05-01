function [Spq] = perfect_shuffle(p,q,for_loop)
% perfect_shuffle produces the a perfect shuffle matrix that turns kron(B,C) into kron(C,B)
%
% ::
%
%
%   [Spq] = perfect_shuffle(p,q)
%   [Spq] = perfect_shuffle(p,q,for_loop)
%
% Args:
%
%    - **p** [numeric]: number of rows of square matrix B
%
%    - **q** [numeric]: number of rows of square matrix C
%
%    - **for_loop** [true|{false}]: if true, a for loop is used. else the
%      procedure is vectorized.
%
% Returns:
%    :
%
%    - **Spq** [matrix]: p*q x p*q matrix such that
%      Spq*kron(B,C)*Spq' = kron(C,B)
%
% Note:
%
% Example:
%
%    m1=3;n1=4; m2=5;n2=7; B=rand(m1,n1);C=rand(m2,n2);BC=kron(B,C);CB=kron(C,B);
%    Sm1_m2=perfect_shuffle(m1,m2); Sn1_n2=perfect_shuffle(n1,n2); max(max(Sm1_m2*BC*Sn1_n2'-CB))==0
%
%    See also:

% References: 
% - Charles F. Van Loan (2000): "The ubiquitous Kronecker product", Journal
%   of Computational and Applied Mathematics 123, pp 85-100.                 
% - Carla Dee Martin (2005): "Higher-order kronecker products and     
%   Tensor decompositions" PhD dissertation, pp 13-14

if nargin<3
    for_loop=false;
end
r=p*q;

Ir=speye(r);

if for_loop
    offset=0;
    Spq=nan(r);
    for irow=1:q
        select_rows=irow:q:r;
        nrows=numel(select_rows);
        Spq(offset+(1:nrows),:)=Ir(select_rows,:);
        offset=offset+nrows;
    end
else
    select_rows=(1:q:r)';
    select_rows=select_rows(:,ones(1,q));
    select_rows=bsxfun(@plus,select_rows,(0:q-1));
    Spq=Ir(select_rows(:),:);
end


end

