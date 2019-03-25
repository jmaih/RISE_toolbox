%  INTERNAL FUNCTION: produces the a perfect shuffle matrix that turns kron(B,C) into kron(C,B)
% 
%  ::
% 
%    [Spq] = perfect_shuffle(p,q)
%    [Spq] = perfect_shuffle(p,q,for_loop)
% 
%  Args:
% 
%     - **p** [numeric]: number of rows of square matrix B
%     - **q** [numeric]: number of rows of square matrix C
%     - **for_loop** [true|{false}]: if true, a for loop is used. else the
%       procedure is vectorized.
% 
%  Returns:
%     :
% 
%     - **Spq** [matrix]: p*q x p*q matrix such that
%       Spq*kron(B,C)*Spq' = kron(C,B)
% 
%  Example:
% 
%     m1=3;n1=4; m2=5;n2=7; B=rand(m1,n1);C=rand(m2,n2);BC=kron(B,C);CB=kron(C,B);
%     Sm1_m2=perfect_shuffle(m1,m2); Sn1_n2=perfect_shuffle(n1,n2); max(max(Sm1_m2*BC*Sn1_n2'-CB))==0
% 
%