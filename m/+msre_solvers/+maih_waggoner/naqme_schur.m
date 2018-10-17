function [X,eigvals,retcode,xtras]=naqme_schur(A,B,C,slvOpts) 
%
% naqme_schur : solves the system A*X^2+B*X+C=0
%
% ::
%
%   [X]=naqme_schur(A,B,C,slvOpts)
%   [X,eigvals,retcode]=naqme_schur(...)
%
% Args:
%
%    - **A** [n x n] : Coefficient matrix on lead variables
%
%     - **B** [n x n] : Coefficient matrix on current variables
%
%     - **C** [n x n] : Coefficient matrix on lagged variables
%
%     - **slvOpts.allSols** [true|{false}] : flag for finding all solutions
%
%     - **slvOpts.checkStab** [true|{false}] : check stability of the system
%
%     - **slvOpts.xplosRoots** [{true}|false] : if all solutions are
%       computed we still can restrict ourselves to solutions that do not
%       involve explosive roots
%
%     - **slvOpts.debug** [true|{false}] : slvOpts.debug or not
%
% Returns:
%    :
%
%     - **X** [n x n x h x k array] : Solution set (k solutions)
%
%     - **eigvals** [empty] : Eigenvalues (Not computed)
%
%     - **retcode** [numeric] : 0 if there is no problem
%
%     - **xtras** [struct] : information on the nature of the different
%       solutions
%
% See also :  udc

% Reference :
% Nicholas J. Higham and Hyun-Min Kim (2000): "Numerical analysis of a
% quadratic matrix equation" IMA Journal of Numerical Analysis (2000) 20,
% 499-519

retcode=0;

xtras=struct();

qz_criterium=sqrt(eps);

if isempty(slvOpts.checkStab),slvOpts.checkStab=false; end

if isempty(slvOpts.xplosRoots),slvOpts.xplosRoots=true; end

if isempty(slvOpts.debug),slvOpts.debug=false; end

if isempty(slvOpts.allSols),slvOpts.allSols=false; end

[mr,m]=size(A);

if mr~=m
    
    error('are you crazy?')
    
end

if slvOpts.debug
    % AX^2+BX+C=0
    test=@(x)A*x^2+B*x+C;
    
end

m=size(A,1);

F=[zeros(m),eye(m)
    -C,-B
    ];

G=[eye(m),zeros(m)
    zeros(m),A];

normalize_by_rows()

[FF, GG, Q, Z] = qz(F,G);

debug1('QZ')

eigvals = ordeig(FF,GG);

% sort in ascending order
%------------------------
[abs_ord_E,order]=sort(abs(eigvals));

% re-order immediately from the most stable to the most unstable
%---------------------------------------------------------------
n=2*m;

clusters=n:-1:1;

klusters(order)=clusters;

[FF,GG,Q,Z] = ordqz(FF,GG,Q,Z,klusters);

eigvals=eigvals(order);%<-- Now same as eigvals = ordeig(FF,GG);

% order
if slvOpts.allSols
    
    if slvOpts.xplosRoots
        
        o=n;
        
    else
        
        o=sum(abs_ord_E<=1+qz_criterium);
        
    end
    
    Sols=set_clusters(o);
    
    xtras.nalt=size(Sols,1);
    
    X=cell(1,xtras.nalt);
    
    xtras.is_stable=false(1,xtras.nalt);
    
    xtras.itersol=0;
    
    for ialt=1:xtras.nalt
        
        clust_i=Sols(ialt,:);
        
        [tmp,stbli,is_bad]=one_solution(clust_i);
        
        ok=~is_bad;
        
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
            
        end
        
    end
    
    X=X(1:xtras.itersol);
    
    xtras.is_stable=xtras.is_stable(1:xtras.itersol);
    
    retcode=zeros(1,xtras.itersol);
    
else
    
    xtras.nalt=1;
    
    xtras.itersol=1;
    
    [X,xtras.is_stable,is_bad]=one_solution(clusters);
        
    if is_bad
        
        retcode=22; % nans in solution
        
    end
    
end

    function [S,is_stable,is_bad]=one_solution(qlusters)
        
        % Re-order by the clusters.
        %--------------------------
        if all(clusters-qlusters==0)
            
            ZZ=Z;
            
        else
            
            [~,~,~,ZZ] = ordqz(FF,GG,Q,Z,qlusters);
            
        end
        
        Z21=ZZ(m+1:end,1:m);
        
        Z11=ZZ(1:m,1:m);
        
        S=Z21/Z11;
        
        is_bad=any(isnan(S(:)))||any(isinf(S(:)));
        
        if is_bad
            
            is_stable=false;
            
        else
            
            if slvOpts.checkStab
                
                is_stable=max(abs(eig(S)))<=1+qz_criterium;
                
            else
                
                is_stable=true;
                
            end
            
        end
        
        if slvOpts.debug
            
            debug2(S)
            
        end
        
    end

    function c=set_clusters(o)
        % here we rotate only the first o out of n elements
        last=clusters(o+1:n);
        
        first=clusters(1:o);
        
        c=cell2mat(utils.gridfuncs.mypermutation(first));
        
        if isempty(last)
            
            return
            
        end
        
        nc=size(c,1);
        
        c=[c,last(ones(nc,1),:)];
        
    end

    function normalize_by_rows()
        
        mynorm=max([max(abs(F),[],2),...
            max(abs(G),[],2)],[],2);
        
        mynorm(mynorm==0)=1;
        
        F=bsxfun(@rdivide,F,mynorm);
        
        G=bsxfun(@rdivide,G,mynorm);
        
    end

    function debug2(sol)
        
        result=max(max(abs(test(sol))))<qz_criterium;
        
        disp(['solution checks :: ',int2str(result)])
        
    end

    function debug1(msg)
        
        if slvOpts.debug
            
            result1=max(max(abs(Q*F*Z-FF)))<qz_criterium;
            
            result2=max(max(abs(Q*G*Z-GG)))<qz_criterium;
            
            result=result1 && result2;
            
            disp([msg,':: ',int2str(result)])
            
        end
        
    end

end