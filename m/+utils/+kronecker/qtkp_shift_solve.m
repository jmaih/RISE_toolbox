function y = qtkp_shift_solve(T,n,c,lambda,alpha)
% qtkp_shift_solve -- Solves the problem (alpha@T1@T2@...@Tp-lambda*I)*y=c
% where @ stands for the kronecker product for quasi-triangular matrices Ti 
%
% ::
%
%	y = qtkp_shift_solve(T,n,c,lambda)
%	y = qtkp_shift_solve(T,n,c,lambda,alpha)
%	y = qtkp_shift_solve(T,[],...)
%
% Args:
%
%    - **T** [cell array]: each element Ti, i=1,2,...,p is a
%      quasi-triangular (square) matrix
%
%    - **n** [vector|{[]}]: number of rows (or columns) for Ti,
%      i=1,2,...,p. 
%
%    - **c** [vector]: right-hand side of the equality
%
%    - **lambda** [scalar]: See formula above
%
%    - **alpha** [scalar|2 x 2 matrix|{1}]: First element of the kronecker
%      product
%
% Returns:
%    :
%
%    - **y** [N x 1 vector]: where N=n1 x n2 x ... x np
%
% See also : utils.kronecker.kp_shift_solve

% References: Carla D. Moravitz Martin and 	Charles F. Van Loan (2006):
% "Shifted Kronecker Product Systems" SIAM Journal on Matrix Analysis and
% Applications archive  Volume 29 Issue 1, December 2006, Pages 184-198

if nargin<5
    
    alpha=1;
    
end

p = numel(T);

if isempty(n)
    
    n=zeros(1,p);
    
    for iii=1:p
        
        n(iii)=size(T{iii},1);
        
    end
    
end

N = prod(n);

if isscalar(alpha)
    
    if alpha~=1
        
        T{1} = alpha*T{1};
        
    end
    
    if p == 1
        
        y = (T{1}-lambda*eye(n(1)))\c;
        
    else
        
        y = zeros(N,1); mp = N/n(1); ii = n(1);
        
        while ii >= 1
            
            if ii > 1 && T{1}(ii,ii-1)~= 0 %(T1 has a 2-by-2 bump)
                
                bump_engine()
                
            else %(T1 does not have a 2-by-2 bump)
                
                no_bump_engine()
                
            end
            
        end
        
    end
    
else % (alpha is 2 x 2)
    
    [Q,S]=schur(alpha,'complex');
    
    nI=size(c,1)/size(Q,1);
        
    d=utils.kronecker.kron_A_I_times_B(Q',c,nI);
    
    z = utils.kronecker.qtkp_shift_solve([{S},T],[2,n(:).'],d,lambda,1);
    
    y=utils.kronecker.kron_A_I_times_B(Q,z,nI);
    
end

    function bump_engine()
        
        idx1=(ii-2)*mp+1:(ii-1)*mp;
        
        idx2=(ii-1)*mp+1:ii*mp;
        
        idx = [idx1,idx2];
        
        y(idx) = utils.kronecker.qtkp_shift_solve(T(2:p),n(2:p),c(idx),...
            lambda,T{1}(ii-1:ii,ii-1:ii));
        
        yy=[y(idx1),y(idx2)];
        
        zz = utils.kronecker.kron_Q1_Qk_times_A(yy,T{2:p});
        
        z1=zz(:,1); z2=zz(:,2);
        
        for jj = 1:ii-2
            
            jdx = (jj-1)*mp+1:jj*mp;
            
            c(jdx) = c(jdx)-T{1}(jj,ii-1)*z1-T{1}(jj,ii)*z2;
            
        end
        
        ii = ii-2;
        
    end

    function no_bump_engine()
        
        idx = (ii-1)*mp+1:ii*mp;
        
        y(idx) = utils.kronecker.qtkp_shift_solve(T(2:p),n(2:p),c(idx),...
            lambda,T{1}(ii,ii));
        
        z = utils.kronecker.kron_Q1_Qk_times_A(y(idx),T{2:p});
        
        for jj = 1:ii-1
            
            jdx = (jj-1)*mp+1:jj*mp;
            
            c(jdx) = c(jdx)-T{1}(jj,ii)*z;
            
        end
        
        ii = ii-1;
        
    end

end