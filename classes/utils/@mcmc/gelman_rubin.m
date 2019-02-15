function [p,R,W,B,Rp]=gelman_rubin(obj,recursive)
% gelman_rubin : computes...
%
% ::
%
%    [p,R,W,B,Rp]=gelman_rubin(obj)
%
%    [p,R,W,B,Rp]=gelman_rubin(obj,recursive)
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

if nargin<2
    
    recursive=false;
    
end

m=obj.nchains;

p=struct();

R=[];

W=[];

B=[];

if m == 1
    
    return
    
end

if recursive
    
    last=2;
    
else
    
    last=obj.npop;
    
end

niter=obj.npop-last+1;

% for istart=last:obj.npop
time_end=obj.start;

R=zeros(obj.nparams,niter); W=R; B=R; V=R;

Rp=zeros(1,niter); aWa=Rp;  aBa=Rp;  aVa=Rp;

% nworkers=utils.parallel.get_number_of_workers();

mm=num2cell(zeros(1,m));

vv=num2cell(zeros(1,m));

batch=obj.draws(1,1).x(:,ones(1,m));
            
for iter=1:niter % parfor(iter=1:niter,nworkers)
    
    for ch=1:m
        
        batch(:,ch)=obj.draws(ch,iter).x;
        
    end
    
    time_end=time_end+1;
    
    [R(:,iter),Wi,Bi,Vi,Rp(iter),aWa(iter),aBa(iter),aVa(iter)]=...
        do_one_recursion();
    
    W(:,iter)=diag(Wi);
    
    B(:,iter)=diag(Bi);
    
    V(:,iter)=diag(Vi);
    
    iter
    
end

for iname=1:obj.nparams
    
    p.(obj.pnames{iname})=struct('within_variance',W(iname,:),...
        'between_variance',B(iname,:),'variance',V(iname,:),'psrf',R(iname,:),...
        'time',obj.start+1:time_end);
    
end

p.multivariate_=struct('within_variance',aWa,...
    'between_variance',aBa,'variance',aVa,'psrf',Rp,...
    'time',obj.start+1:time_end);


    function [R,W,B,V,Rp,aWa,aBa,aVa]=do_one_recursion()
                        
        n=iter;
                
        % 2. Discard half of the draws
        %-----------------------------
%         newDraws=load_draws();
        
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
                
                w=w+im*newcov();
                
            end
            
            function cc=newcov()
                
                d=batch(:,j);
                
                [mm{j},vv{j},n]=utils.moments.recursive(mm{j},vv{j},d,iter);
                
                cc=vv{j};
                
            end
            
        end
        
        function b=between_chain_variance()
            
            theta_bar=zeros(obj.nparams,m);
            % use the results already computed in the within chain
            for jj=1:m
                
                theta_bar(:,jj)=mm{jj};
                
            end
            
            b=n*cov(theta_bar.');
                                    
        end
                
    end

end
