%  INTERNAL FUNCTION: Produces the a perfect shuffle matrix that turns kron(B,C) into kron(C,B)
% 
%  ::
% 
%    [CB] = shuffle(BC,bdim,cdim)
% 
%  Args:
% 
%     - **BC** [matrix]: kron(B,C) of size bdim(1)*cdim(1) x bdim(2)*cdim(2)
%     - **bdim** [scalar|2-element vector]: dimensions of matrix B. Could be a
%       scalar if B is square
%     - **cdim** [scalar|2-element vector]: dimensions of matrix C. Could be a
%       scalar if C is square
% 
%  Returns:
%     :
% 
%     - **CB** [matrix]: kron(C,B) bdim(1)*cdim(1) x bdim(2)*cdim(2) matrix
% 
%  Note:
% 
%     In some applications, we know kron(B,C) but cannot compute kron(C,B)
%     directly or it can be expensive to compute.
% 
%  Example:
% 
%     ::
% 
%        m1=3;n1=4; m2=5;n2=7; B=rand(m1,n1);C=rand(m2,n2);BC=kron(B,C);
%        CB=shuffle(BC,[m1,n1],[m2,n2]); max(max(kron(C,B)-CB))==0
% 
%        m1=3; m2=5; B=rand(m1); C=rand(m2); BC=kron(B,C);
%        CB=shuffle(BC,m1,m2); max(max(kron(C,B)-CB))==0
% 
%