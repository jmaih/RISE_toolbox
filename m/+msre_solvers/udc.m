function [X,eigvals,retcode,xtras]=udc(A,B,C,allsols,msv_only,debug)
%
% udc Undetermined coefficients solution algorithm for DSGE models.
% The procedure can find all possible solutions for a constant-parameter
% DSGE model or a regime-switching DSGE model with diagonal transition
% matrix.
%
% ::
%
%   [Tz_pb,eigvals,retcode]=udc(A,B,C,Q,T0,TolFun,maxiter)
%
% Args:
%
%    - **A** [n x n] : Coefficient matrix on lead variables
%
%     - **B** [n x n] : Coefficient matrix on current variables
%
%     - **C** [n x n] : Coefficient matrix on lagged variables
%
%     - **allsols** [true|{false}] : flag for finding all solutions
%
%     - **msv_only** [{true}|false] : return only MSV solutions
%
%     - **debug** [true|{false}] : debug or not
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
% See also :  groebner

if nargin<6
    
    debug=[];
    
    if nargin<5
        
        msv_only=[];
        
        if nargin<4
            
            allsols=[];
            
        end
        
    end
    
end

retcode=0;

xtras=struct();

qz_criterium=sqrt(eps);

if isempty(msv_only),msv_only=true; end

if isempty(debug),debug=false; end

if isempty(allsols),allsols=false; end

[mr,m]=size(A);

if mr~=m
    
    error('are you crazy?')
    
end

msv_cols=all(C==0);

XI=[-B,-C
    eye(m),zeros(m)];

DELTA=[A,zeros(m)
    zeros(m),eye(m)];

[V,D]=eig(XI,DELTA);

debug1('abs(XI*V-DELTA*V*D)<1e-12')

eigvals=diag(D);

[~,order]=sort(abs(eigvals));

V=V(:,order);

D=D(order,order);

debug1('sorted V and D:: abs(XI*V-DELTA*V*D)<1e-12')

V2=V(m+1:end,:);

debug2()

% avoid computations with infinite eigenvalues
%---------------------------------------------
bad=~isfinite(diag(D));

V2=V2(:,~bad);

D=D(~bad,~bad);

eigvals=diag(D);

n=size(D,1);

if n<m
    
    error('No solution: I don''t expect this case to ever happen')
    
end

if allsols
    
    Sols=utils.gridfuncs.ordered_combos_without_repetition(1:n,m);
    
    xtras.nalt=size(Sols,1);
    
    X=cell(1,xtras.nalt);
    
    xtras.is_stable=false(1,xtras.nalt);
    
    xtras.is_msv=false(1,xtras.nalt);
    
    xtras.itersol=0;
    
    for ialt=1:xtras.nalt
        
        srange=Sols(ialt,:);
        
        [tmp,stbli,is_msvi,is_bad]=one_solution(srange);
        
        ok=(~msv_only||(msv_only && is_msvi))&&~is_bad;
        
        if ok
            
            for ialt2=1:xtras.itersol-1
                
                ok=max(abs(tmp(:)-X{ialt2}(:)))>1e-5;
                
                if ~ok
                    
                    break
                    
                end
                
            end
            
        end
        
        if ok
            
            xtras.itersol=xtras.itersol+1;
            
            X{xtras.itersol}=tmp;
            
            xtras.is_stable(xtras.itersol)=stbli;
            
            xtras.is_msv(xtras.itersol)=is_msvi;
            
        end
        
    end
    
    X=X(1:xtras.itersol);
    
    xtras.is_stable=xtras.is_stable(1:xtras.itersol);
    
    xtras.is_msv=xtras.is_msv(1:xtras.itersol);
    
else
    
    xtras.nalt=1;
    
    xtras.itersol=1;
    
    [X,xtras.is_stable,xtras.is_msv,is_bad]=one_solution(1:m);
    
    X=real(X);
    
    if is_bad
        
        retcode=22; % nans in solution
        
    end
    
end

    function [S,is_stable,is_msv,is_bad]=one_solution(range)
        
        OMG=V2(:,range);
        
        LAMB=D(range,range);
        
        S=OMG*LAMB/OMG;
        
        is_bad=any(isnan(S(:)))||any(isinf(S(:)));
        
        if is_bad
            
            is_msv=false; is_stable=false;
            
        else
            
            is_msv=all(all(S(:,msv_cols)==0)); % I don't know whether this is too strict
            
            is_stable=max(abs(diag(LAMB)))<=1+qz_criterium;
            
        end
        
        if debug
            
            result=max(max(abs(A*S^2+B*S+C)))<qz_criterium;
            
            disp(['A*X^2+B*X+C=0:: ',int2str(result)])
            
            if ~result
                
                keyboard
                
            end
            
        end
        
    end

    function debug1(msg)
        
        if debug
            
            result=max(max(abs(XI*V-DELTA*V*D)))<1e-12;
            
            disp([msg,':: ',int2str(result)])
            
        end
        
    end

    function debug2()
        
        if debug
            
            V1=V(1:m,:);
            
            tmp=V1./V2;
            
            result=max(abs(tmp(1,:)-diag(D).'))<1e-12;
            
            disp(['V2=lamb*V1:: ',int2str(result)])
            
        end
        
    end

end