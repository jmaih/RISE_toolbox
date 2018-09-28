function x=kp_shift_solve(A,lambda,b,do_real)
% kp_shift_solve -- Solves the problem (A1@A2@...@Ap-lambda*I)*x=b
% where @ stands for the kronecker product
%
% ::
%
%	x = kp_shift_solve(A,lambda,b)
%	x = kp_shift_solve(A,lambda,b,do_real)
%
% Args:
%
%    - **A** [cell array]: each element Ai, i=1,2,...,p is a
%      square matrix with number of rows ni
%
%    - **b** [vector]: right-hand side of the equality
%
%    - **lambda** [scalar]: See formula above
%
%    - **do_real** [true|{false}]: See formula above
%
% Returns:
%    :
%
%    - **x** [N x 1 vector]: where N=n1 x n2 x ... x np
%
% See also : utils.kronecker.qtkp_shift_solve

% References: Carla D. Moravitz Martin and 	Charles F. Van Loan (2006):
% "Shifted Kronecker Product Systems" SIAM Journal on Matrix Analysis and
% Applications archive  Volume 29 Issue 1, December 2006, Pages 184-198

if nargin<4
    
    do_real=[];
    
end

if isempty(do_real)
    
    do_real=false;
    
end

if do_real
    
    type='real';
    
else
    
    type='complex';
    
end

p=numel(A);

n=zeros(1,p);

T=cell(1,p);

U=T;

Uprime=U;

for ii=1:p
    
    n(ii)=size(A{ii},1);
    
    [U{ii},T{ii}] = schur(A{ii},type);
    
    Uprime{ii}=U{ii}';
    
end

c=utils.kronecker.kron_Q1_Qk_times_A(b,Uprime{:});

y = utils.kronecker.qtkp_shift_solve(T,n,c,lambda);

x=utils.kronecker.kron_Q1_Qk_times_A(y,U{:});

if ~do_real
    
    x=real(x);
    
end

end
