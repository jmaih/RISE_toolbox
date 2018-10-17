function [Tz_pb,eigvals,retcode]=dsge_schur(Alead,Acurr,Alag,Q,~,TolFun,...
    maxiter,varargin)
%
% dsge_schur : Schur solution algorithm for DSGE models.
% The procedure can find all possible solutions for a constant-parameter
% DSGE model or a regime-switching DSGE model with diagonal transition
% matrix.
%
% ::
%
%   [Tz_pb,eigvals,retcode]=dsge_schur(Alead,Acurr,Alag,Q,T0,TolFun,maxiter)
%
%   [Tz_pb,eigvals,retcode]=dsge_schur(Alead,Acurr,Alag,Q,T0,TolFun,maxiter,varargin)
%
% Args:
%
%    - **Alead** [n x n x h x h array] : The suffix 01 is there to
%      indicate that the Aplus matrices are multiplied with the transition
%      probabilities: Jacobian of lead variables
%
%     - **Acurr** [n x n x h array] : Jacobian of contemporaneous variables in
%       each regime
%
%     - **Alag** [n x n x h array] : Jacobian of lagged variables in each
%       regime
%
%     - **Q** [h x h matrix] : Transition matrix (used only if refinement)
%
%     - **T0** [n x n x h array] : Initial guess for the solution (Not used)
%
%     - **TolFun** [numeric] : Tolerance criterion for solution used in the
%       refinement of the solution
%
%     - **maxiter** [numeric] : Maximum number of iterations used in the
%       refinement of the solution
%
%     - **varargin** [] : additional arguments
%        - **refine** [true|{false}] : solve for sigma=1 instead of just
%          sigma=0
%        - **checkStab** [true|{false}] : check the stability of the
%           solution by the eigenvalues
%        - **allSols** [true|{false}] : flag for finding all solutions
%        - **msvOnly** [{true}|false] : return only the minimum state
%          variable solutions
%        - **xplosRoots** [{true}|false] : if all solutions are
%          computed we still can restrict ourselves to solutions that do not
%          involve explosive roots
%        - **debug** [true|{false}] : debug or not
%
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
% See also :  dsge_groebner, dsge_udc


slvOpts=msre_solvers.maih_waggoner.set_solve_options(varargin{:});

[n,h]=reformat_inputs();

Alead2=zeros(n*h);

Acurr2=zeros(n*h);

Alag2=zeros(n*h);

for s0=1:h
    
    row_range=(s0-1)*n+(1:n);
    
    Acurr2(row_range,row_range)=Acurr{s0};
    
    Alag2(row_range,row_range)=Alag{s0};
    
    for s1=1:h
        
        As0s1=Alead{s0,s1};
        
        if s0==s1
            
            Alead2(row_range,row_range)=As0s1;
            
        else
            
            if max(abs(As0s1(:)))>1e-10
                
                Tz_pb=[];
                
                eigvals=[];
                
                retcode=28; % error('transition matrix must be diagonal')
                
                return
                
            end
            
        end
        
    end
    
end

if slvOpts.allSols
    
    doall=true;
    
    [X,eigvals,retcode]=msre_solvers.maih_waggoner.naqme_schur(Alead2,...
        Acurr2,Alag2,slvOpts);
    
else
    
doall=false;

    eigvals=cell(h,1);
    
    X=zeros(n*h);
    
    retcode=zeros(1,h);
    
    for s0=1:h
        
        [Xs0,eigvals{s0},retcode(s0)]=msre_solvers.maih_waggoner.naqme_schur(...
            full(Alead{s0,s0}),full(Acurr{s0}),full(Alag{s0}),slvOpts);
        
        if retcode(s0)
            
            break
            
        end
        
        row_range=(s0-1)*n+(1:n);
        
        X(row_range,row_range)=Xs0;
        
    end
    
    retcode=max(retcode);
    
    if all(retcode)
        
        eigvals=[];
        
        Tz_pb=[];
        
        return
        
    end
        
end

if ~iscell(X)
    
    X={X};
    
end

k=numel(X);

Tz_pb=zeros(n,n,h,k);

Sproto=zeros(n,n,h);

refine_data=[];

if slvOpts.refine
    
    opts=msre_solvers.maih_waggoner.update_options(TolFun,maxiter,slvOpts.debug);

    Alead_diag=msre_solvers.maih_waggoner.diag_cell(Alead);
    
    refine_data={Acurr,Alead_diag,Q,opts};
    
end

for isol=1:k
    
    if ~retcode(isol)
        
        [Tz_pb(:,:,:,isol),retcode(isol)]=...
            msre_solvers.maih_waggoner.set_one_solution(X{isol},Sproto,...
            doall,refine_data);
        
    end
    
end

    function [n,h]=reformat_inputs()
        
        [n,nc,h,h2]=size(Alead);
        
        if n~=nc, error('wrong size of first input'), end
        
        if h~=h2, error('wrong size of first input'), end
        
        tmp_lead=Alead; tmp_curr=Acurr; tmp_lag=Alag;
        
        Alead=cell(h); Acurr=cell(1,h); Alag=cell(1,h);
        
        for r0=1:h
            
            Acurr{r0}=sparse(tmp_curr(:,:,r0));
            
            Alag{r0}=sparse(tmp_lag(:,:,r0));
            
            for r1=1:h
                
                Alead{r0,r1}=sparse(tmp_lead(:,:,r0,r1));
                
            end
            
        end
        
    end

end
