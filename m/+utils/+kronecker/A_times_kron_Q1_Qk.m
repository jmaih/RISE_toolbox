function C=A_times_kron_Q1_Qk(A,varargin)
% A_times_kron_Q1_Qk -- computes A(Q1*Q2*...*Qk) where * denotes the
% kronecker product
%
% ::
%
%
% - C=A_times_kron_Q1_Qk(A,Q1,Q2,...,Qk)
%
% Args:
%
%    - **A** [matrix]: matrix on the left-hand side
%
%    - **Qi** [matrix]: matrix on the kronecker block
%
% Returns:
%    :
%
%    - **C** [matrix]: result
%
% Note:
%
% Example:
%
%    kronall=@utils.kronecker.kronall;
%    rb=3;cb=4;B=rand(rb,cb);
%    rc=2;cc=7;C=rand(rc,cc);
%    rd=4;cd=6;D=rand(rd,cd);
%    re=2;ce=5;E=rand(re,ce);
%    ra=3;ca=rb*rc*rd*re;A=rand(ra,ca);
%    H=A*kronall(B,C,D,E);
%    Htest=utils.kronecker.A_times_kron_Q1_Qk(A,B,C,D,E);
%    max(abs(Htest(:)-H(:)))
%
%    See also:

k=length(varargin);

rq=zeros(1,k);

cq=rq;

for iq=1:k
    
    [rq(iq),cq(iq)]=size(varargin{iq});

end

[~,ca]=size(A);

for iq=1:k
    
    Qi=varargin{iq};
    
    if iq==1
        
        r0=ca/rq(1);
        
        C=utils.kronecker.A_times_kron_B_I(A,Qi,r0);
    
    elseif iq==k
        
        rk=prod(cq(1:k-1));
        
        C=utils.kronecker.A_times_kron_I_B(C,Qi,rk);
    
    else
        
        r_right=prod(rq(iq+1:end));
        
        r_left=size(C,2)/(rq(iq)*r_right);
        
        C=utils.kronecker.A_times_kron_I_B_I(C,Qi,r_left,r_right);
    
    end
    
end