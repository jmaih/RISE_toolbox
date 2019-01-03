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

my_recursion=@do_one_recursion;

R=zeros(obj.nparams,niter); W=R; B=R; V=R;

Rp=zeros(1,niter); aWa=Rp;  aBa=Rp;  aVa=Rp;

nworkers=utils.parallel.get_number_of_workers();

parfor(iter=1:niter,nworkers)
    
    time_end=time_end+1;
    
    [R(:,iter),W(:,iter),B(:,iter),V(:,iter),Rp(iter),aWa(iter),aBa(iter),aVa(iter)]=...
        my_recursion(obj.draws(:,1:last+iter-1)); %#ok<PFBNS>
    
end

for iname=1:obj.nparams
    
    p.(obj.pnames{iname})=struct('within_variance',W(iname,:),...
        'between_variance',B(iname,:),'variance',V(iname,:),'psrf',R(iname,:),...
        'time',obj.start+1:time_end);
    
end

p.multivariate_=struct('within_variance',aWa,...
    'between_variance',aBa,'variance',aVa,'psrf',Rp,...
    'time',obj.start+1:time_end);


    function [R,W,B,SIG2,Rp,aWa,aBa,aVa]=do_one_recursion(these_draws)
        
        two_n=size(these_draws,2);
        
        %         if obj.i_dropped>0
        %             % do not drop more than necessary
        %             n=two_n;
        %
        %         else
        
        n=floor(0.5*two_n);
        
        %         end
        
        % 2. Discard half of the draws
        %-----------------------------
        newDraws=load_draws();
        
        % 3. Calculate the within-chain and between-chain variances
        %----------------------------------------------------------
        W=within_chain_variance();
        
        B=between_chain_variance();
        
        % 4. calculate the estimated variance of the parameter as a weighted sum of
        % the within-chain and between-chain variances
        SIG2=(1-1/n)*W+1/n*B;
        
        %         Vhat=SIG2+B/(m*n);
        
        %         R=(SIG2+B/(m*n))/W
        
        % 5. Calculate the potential scale reduction factor
        %------------------------------ --------------------
        R=sqrt(SIG2./W);
        
        % 5. Calculate the Multivariate potential scale reduction factor
        %------------------------------ --------------------------------
        
        [eigVect,D]=eig((W\B)/n);
        
        L=diag(D);
        
        L1=max(L);
        
        best=L==L1;
        
        a=eigVect(:,best);
        
        V_=(1-1/n)*W+(1+1/m)*1/n*B;
        
        aVa=a.'*V_*a;
        
        aWa=a.'*W*a;
        
        aBa=a.'*B*a;
        
        Rp=(n-1)/n+(m+1)/m*L1;
        
        function [w,s2j]=within_chain_variance()
            
            s2j=zeros(obj.nparams,m);
            
            for ichain=1:m
                
                s2j(:,ichain)=var(newDraws{ichain},1,2);
                
            end
            
            w=mean(s2j,2);
            
        end
        
        function b=between_chain_variance()
            
            theta_bar=zeros(obj.nparams,m);
            
            for ichain=1:m
                
                theta_bar(:,ichain)=mean(newDraws{ichain},2);
                
            end
            
            theta_bar_bar=mean(theta_bar,2);
            
            b=n/(m-1)*sum(bsxfun(@minus,theta_bar,theta_bar_bar).^2,2);
            
        end
        
        function newDraws=load_draws()
            
            newDraws=these_draws(:,end-n+1:end);
            
            tmp=cell(m,1);
            
            for ichain=1:m
                
                tmp{ichain}=[newDraws(ichain,:).x];
                
            end
            
            newDraws=tmp;
            
        end
        
    end

end
