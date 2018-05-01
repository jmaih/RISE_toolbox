function [CB] = shuffle(BC,bdim,cdim)
% perfect_shuffle produces the a perfect shuffle matrix that turns kron(B,C) into kron(C,B)
%
% ::
%
%
%   [CB] = shuffle(BC,bdim,cdim)
%
% Args:
%
%    - **BC** [matrix]: kron(B,C) of size bdim(1)*cdim(1) x bdim(2)*cdim(2)
%
%    - **bdim** [scalar|2-element vector]: dimensions of matrix B. Could be a
%      scalar if B is square
%
%    - **cdim** [scalar|2-element vector]: dimensions of matrix C. Could be a
%      scalar if C is square
%
% Returns:
%    :
%
%    - **CB** [matrix]: kron(C,B) bdim(1)*cdim(1) x bdim(2)*cdim(2) matrix
%
% Note:
%
%    In some applications, we know kron(B,C) but cannot compute kron(C,B)
%    directly or it can be expensive to compute.
%
% Example:
%
%    m1=3;n1=4; m2=5;n2=7; B=rand(m1,n1);C=rand(m2,n2);BC=kron(B,C);
%    CB=shuffle(BC,[m1,n1],[m2,n2]); max(max(kron(C,B)-CB))==0
%
%    m1=3; m2=5; B=rand(m1); C=rand(m2); BC=kron(B,C);
%    CB=shuffle(BC,m1,m2); max(max(kron(C,B)-CB))==0
%
%    See also:

% References: 
% - Charles F. Van Loan (2000): "The ubiquitous Kronecker product", Journal
%   of Computational and Applied Mathematics 123, pp 85-100.                 
% - Carla Dee Martin (2005): "Higher-order kronecker products and     
%   Tensor decompositions" PhD dissertation, pp 13-14

if numel(bdim)==1
    bdim=bdim*ones(1,2);
end

if numel(cdim)==1
    cdim=cdim*ones(1,2);
end

m1=bdim(1);
n1=bdim(2);
m2=cdim(1);
n2=cdim(2);

% BC=sparse(BC);

[mnr,mnc]=size(BC);

if mnr~=(m1*m2)||mnc~=(n1*n2)
    error('matrix size does not match the inputed dimensions')
end

Sm1_m2=utils.kronecker.perfect_shuffle(m1,m2);

CB=Sm1_m2*BC;
if m1==n1 && m2==n2
    CB=CB*Sm1_m2';
else
    Sn1_n2=utils.kronecker.perfect_shuffle(n1,n2);
    CB=CB*Sn1_n2';
end

end

