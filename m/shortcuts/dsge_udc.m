function [Tz_pb,eigvals,retcode]=dsge_udc(Alead,Acurr,Alag,Q,~,TolFun,maxiter,varargin)
%
% dsge_udc undetermined coefficients-type of solution algorithm for DSGE models.
% The procedure can find all possible solutions for a constant-parameter
% DSGE model or a regime-switching DSGE model with diagonal transition
% matrix.
%
% ::
%
%   [Tz_pb,eigvals,retcode]=dsge_udc(Alead,Acurr,Alag,Q,T0,TolFun,maxiter)
%
%   [Tz_pb,eigvals,retcode]=dsge_udc(Alead,Acurr,Alag,Q,T0,TolFun,maxiter,varargin)
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
%     - **Q** [h x h matrix] : Transition matrix (Not used)
%
%     - **T0** [n x n x h array] : Initial guess for the solution (Not used)
%
%     - **TolFun** [numeric] : Tolerance criterion for solution (Not used)
%
%     - **maxiter** [numeric] : Maximum number of iterations (Not used)
%
%     - **varargin** [] : additional arguments
%        - **allsols** [true|{false}] : flag for finding all solutions
%        - **msv_only** [{true}|false] : return only MSV solutions
%        - **debug** [true|{false}] : debug or not
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
% See also :  dsge_groebner

refine=false;

debug=false;

larg=length(varargin);

if larg>3
    
    debug=varargin{3};
    
end

opts=struct();

update_options()

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


doall=false;

if isempty(varargin)||~varargin{1} % not all solutions
    
    eigvals=cell(h,1);
    
    X=zeros(n*h);
    
    retcode=zeros(1,h);
    
    for s0=1:h
        
        [Xs0,eigvals{s0},retcode(s0)]=msre_solvers.udc(...
            full(Alead{s0,s0}),full(Acurr{s0}),full(Alag{s0}),varargin{:});
        
        if retcode(s0)
            
            break
            
        end
        
        row_range=(s0-1)*n+(1:n);
        
        X(row_range,row_range)=Xs0;
        
    end
    
    retcode=max(retcode);
    
    if retcode
        
        eigvals=[];
        
        Tz_pb=[];
        
        return
        
    end
    
else
    
    doall=true;
    
    [X,eigvals,retcode]=msre_solvers.udc(Alead2,Acurr2,Alag2,varargin{:});
    
end

if ~iscell(X)
    
    X={X};
    
end

k=numel(X);

Tz_pb=zeros(n,n,h,k);

Sproto=zeros(n,n,h);

if refine
    
    Alead_diag=diag_cell(Alead);
    
end

for isol=1:k
    
    Tz_pb(:,:,:,isol)=set_one_solution(X{isol});
    
end

    function [S,itercode,retcode]=set_one_solution(X)
        
        itercode=0;
        
        [X,retcode]=dsge_realize_solution(X,doall);
        
        S=Sproto;
        
        for ireg=1:h
            
            row_range=(ireg-1)*n+(1:n);
            
            S(:,:,ireg)=X(row_range,row_range);
            
        end
        
        if refine
            
            [T2,itercode,retcode]=solve_pretzel(S,Acurr,Alead_diag,Q,opts);
            
            S=S+T2;
            
        end
        
    end

    function update_options()
        
        opts.fix_point_TolFun=TolFun;
        
        opts.fix_point_maxiter=maxiter;
        
        opts.fix_point_verbose=debug;
        
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

function Cout=diag_cell(Cin)

[h,h2]=size(Cin);

if h~=h2
    
    error('cell array should be square')
    
end

Cout=cell(1,h);

for ireg=1:h
    
    Cout{ireg}=Cin{ireg,ireg};
    
end

end