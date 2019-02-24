function [p,R,W,B,Rp]=gelman_rubin(obj,~)
% gelman_rubin : computes...
%
% ::
%
%    [p,R,W,B,Rp]=gelman_rubin(obj)
%
% Args:
%
%    obj (mcmc object): mcmc object
%
%    recursive (true|{false}):
%
% Returns:
%    :
%
%    - **p** (struct): with parameter names as fields and for each
%    field a structure containing information on PSRF, Within and Between
%    variances and the variance. N.B: One of the fields name is
%    "multivariate_" and it represents the aggregated statistics.
%
%    - **R** (Matrix): Potential scale reduction factor for each parameter
%
%    - **W** (Matrix): Within variance
%
%    - **B** (matrix): Between variance
%
%    - **Rp** (vector): Multivariate Potential Scale Reduction Factor
%
% Warning:
%
%    - This function requires multiple chains of MCMC samples. See
%      **nchain** option of samplers.
%
% References:
%
%    - :cite:`gelman1992inference`
%

m=obj.nchains;

p=struct();

R=[];

W=[];

B=[];

if m == 1
    
    return
    
end

Origin=1000;

ndraws=obj.npop;

% So that the computational time does not increase too much with the number of simulations.
StepLength = 100; % StepLength = ceil((ndraws-Origin)/100);

recurGrid = (Origin:StepLength:ndraws);

if isempty(recurGrid)
    
    recurGrid = ndraws;
    
elseif recurGrid(length(recurGrid))<ndraws
    
    recurGrid = [recurGrid ndraws];
    
end

nsteps=numel(recurGrid);

R=zeros(obj.nparams,nsteps); W=R; B=R; V=R;

Rp=zeros(1,nsteps); aWa=Rp;  aBa=Rp;  aVa=Rp;

% place holders for moments
%---------------------------
mm=num2cell(zeros(1,m));

vv=mm;

prev_n=0;

for istep=1:nsteps 
        
    n=recurGrid(istep);
    
    subdraws=obj.draws(:,prev_n+1:n);
    
    [R(:,istep),Wi,Bi,Vi,Rp(istep),aWa(istep),aBa(istep),aVa(istep)]=...
        do_one_recursion();
    
    W(:,istep)=diag(Wi);
    
    B(:,istep)=diag(Bi);
    
    V(:,istep)=diag(Vi);
    
    prev_n=n;
    
end

for iname=1:obj.nparams
    
    p.(obj.pnames{iname})=struct('within_variance',W(iname,:),...
        'between_variance',B(iname,:),'variance',V(iname,:),'psrf',R(iname,:),...
        'time',1:nsteps);
    
end

p.multivariate_=struct('within_variance',aWa,...
    'between_variance',aBa,'variance',aVa,'psrf',Rp,...
    'time',1:nsteps);

    function [R,W,B,V,Rp,aWa,aBa,aVa]=do_one_recursion()
                        
        newDraws=load_draws();
        
        % 3. Calculate the within-chain and between-chain variances
        %----------------------------------------------------------
        W=within_chain_variance();
        
        B=between_chain_variance();
        
        % 4. calculate the estimated variance of the parameter as a weighted sum of
        % the within-chain and between-chain variances
        V=(n-1)/n*W+(1+1/m)*B/n;
                
        % 5. Calculate the potential scale reduction factor
        %------------------------------ --------------------
        R=sqrt(diag(V./W));
        
        % 5. Calculate the Multivariate potential scale reduction factor
        %------------------------------ --------------------------------
        
        iWB=W\B;
        
        if any(isnan(iWB(:)))||~all(isfinite(iWB(:)))
            
            aVa=nan; aWa=nan; aBa=nan; Rp=nan;
            
            return
            
        end
        
        [eigVect,D]=eig(iWB/n);
        
        L=diag(D);
        
        L1=max(L);
        
        best=L==L1;
        
        a=eigVect(:,best);
                
        aVa=a.'*V*a;
        
        aWa=a.'*W*a;
        
        aBa=a.'*B*a;
        
        Rp=(n-1)/n+(m+1)/m*L1;
        
        function [w]=within_chain_variance()
            
            w=0;
            
            im=1/m;
            
            for j=1:m
                
                update_moments()
                
                w=w+im*vv{j};
                
            end
            
            
            function update_moments()
                
                iter=prev_n+1;
                
                [mm{j},vv{j}]=utils.moments.recursive(mm{j},vv{j},newDraws{j},iter,prev_n);
                
                % verify with 
                % mean([obj.draws(j,1:n).x],2)-mm{j}
                % cov([obj.draws(j,1:n).x].',1)-vv{j}

            end
            
        end
        
        function b=between_chain_variance()
            
            theta_bar=zeros(obj.nparams,m);
            
            for jj=1:m
                
                theta_bar(:,jj)=mm{jj};
                
            end
            
            b=n*cov(theta_bar.');
                                    
        end
        
        function newDraws=load_draws()
                        
            tmp=cell(m,1);
            
            for ichain=1:m
                
                tmp{ichain}=[subdraws(ichain,:).x];
                
            end
            
            newDraws=tmp;
            
        end
        
    end

end
