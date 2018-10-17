function [Tz_pb,eigvals,retcode]=groebner(Gplus01_0,A0_0,A_0,~,~,~,~,...
    minimal)
% 
% groebner Groebner basis-based solution for DSGE models. Symbolic
% solution of algebraic equations implied by a DSGE model. The procedure
% attempts to find all possible solutions.
%
% ::
%
%   [Tz_pb,eigvals,retcode]=groebner(Gplus01_0,A0_0,A_0,Q,T0,TolFun,maxiter)  
%
% Args:
%
%    - **Gplus01_0** [n x n x h x h array] : The suffix 01 is there to
%      indicate that the Aplus matrices are multiplied with the transition
%      probabilities: Jacobian of lead variables
%
%     - **A0_0** [n x n x h array] : Jacobian of contemporaneous variables in
%       each regime
%
%     - **A_0** [n x n x h array] : Jacobian of lagged variables in each
%       regime 
%
%     - **Q** [h x h matrix] : Transition matrix (Not used) 
%
%     - **T0** [n x n x h array] : Initial guess for the solution (Not used) 
%
%     - **TolFun** [numeric] : Tolerance criterion for solution (Not used) 
%
%     - **maxiter** [numeric] : Maximum number of iterations (Not used) 
%
% Returns:
%    :
%
%     - **Tz_pb** [n x n x h x k array] : Solution set (k solutions)
%
%     - **eigvals** [empty] : Eigenvalues (Not computed) 
%
%     - **retcode** [numeric] : 0 if there is no problem 
%
% N:B:
%    :
%
%    This solver is not guaranteed to give solutions and its computational
%    costs increase rapidly with the size of the problem. For instance I
%    have not been able to solve a problem with 18 unknowns. It is
%    therefore advisable to have a problem as small as possible.
%
%    If Matlab shows busy for a long time, it is better just to kill it
%    rather than waiting.
%
%    This procedures requires Matlab version 2018a or above.
%
% See also :

if nargin<8
    
    minimal=[];
    
end

if isempty(minimal),minimal=true; end

if minimal
    
    aminus=any(any(A_0,1),3);
    
    aplus=any(any(any(Gplus01_0,1),3),4);
    
end

eigvals=[];

h=size(A0_0,3);

Ap=cell(h);

A0=cell(1,h);

A_=cell(1,h);

for s0=1:h
    
    for s1=1:h
        
        if minimal
            
            Ap{s0,s1}=Gplus01_0(:,aplus,s0,s1);
            
        else
            
            Ap{s0,s1}=Gplus01_0(:,:,s0,s1);
            
        end
        
    end
    
    A0{s0}=A0_0(:,:,s0);
    
    if minimal
        
        A_{s0}=A_0(:,aminus,s0);
        
    else
        
        A_{s0}=A_0(:,:,s0);
        
    end
    
end

if minimal
    
    T=utils.solve.grobner_solve(Ap,A0,A_,aplus,aminus);
    
else
    
    T=utils.solve.grobner_solve(Ap,A0,A_);
    
end

n=numel(aplus);

nsol=numel(T);

Tz_pb=zeros(n,n,h,nsol);

retcode=zeros(1,nsol);

for isol=1:nsol
    
    if minimal
        
        Tz_pb(:,aminus,:,isol)=T{isol};
        
    else
        
        Tz_pb(:,:,:,isol)=T{isol};
        
    end
    
    if ~utils.error.validComplex(T{isol})% any(isnan(T{isol}(:)))||any(isinf(T{isol}(:)))
        
        retcode(isol)=22;
        
    end
    
end

end