%  INTERNAL FUNCTION: computes A(Q1*Q2*...*Qk) where * denotes the kronecker product
% 
%  ::
% 
%     C=A_times_kron_Q1_Qk(A,Q1,Q2,...,Qk)
% 
%  Args:
% 
%     - **A** [matrix]: matrix on the left-hand side
%     - **Qi** [matrix]: matrix on the kronecker block
% 
%  Returns:
%     :
% 
%     - **C** [matrix]: result
% 
%  Example:
% 
%     ::
% 
%        kronall=@utils.kronecker.kronall;
%        rb=3;cb=4;B=rand(rb,cb);
%        rc=2;cc=7;C=rand(rc,cc);
%        rd=4;cd=6;D=rand(rd,cd);
%        re=2;ce=5;E=rand(re,ce);
%        ra=3;ca=rb*rc*rd*re;A=rand(ra,ca);
%        H=A*kronall(B,C,D,E);
%        Htest=utils.kronecker.A_times_kron_Q1_Qk(A,B,C,D,E);
%        max(abs(Htest(:)-H(:)))
% 
%