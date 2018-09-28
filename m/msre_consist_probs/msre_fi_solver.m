function [T,retcode,c,Q]=msre_fi_solver(obj,inner_probs,T0,x0,control_shocks,options)

if nargin < 6
    
    options=struct();
    
    if nargin<5
        
        control_shocks = [];
        
        if nargin < 4
            
            x0 = []; % ordered according to the order of m.endogenous.name
            
            if nargin < 3
                
                T0 = [];
                
            end
            
        end
        
    end
    
end

endo_nbr=obj.endogenous.number;

exo_nbr=sum(obj.exogenous.number);

regime_nbr=obj.markov_chains.regimes_number;

if isempty(T0);
    
    T0.Tx=zeros(endo_nbr,endo_nbr,regime_nbr);
    
    T0.Tsig=zeros(endo_nbr,1,regime_nbr);
    
    T0.Te=zeros(endo_nbr,exo_nbr,regime_nbr);
    
    T0.ss=[];
    
end

if isempty(control_shocks)
    
    control_shocks=randn(exo_nbr,1,10000);
    
end

order_var=obj.order_var;

inner_probs=subutils.process_probabilities(obj,inner_probs);

aie=speye(endo_nbr);

% Place holders for constant state matrices
own_stuff=[];

nxx=numel(T0.Tx);

[obj,structural_matrices,retcode]=compute_steady_state(obj);

ss=obj.solution.ss;

for ireg=1:numel(ss)
    
    ss{ireg}=ss{ireg}(order_var);
    
end

x0=x0(order_var);

[t,~,retcode]=fix_point_iterator(@iterator,subutils.vectorizer(T0),...
    options);

T=struct();

[T.Tx,T.Tsig,T.Te,T.ss]=resizer(t);

    function [Tx,Tsig,Te,ss1]=resizer(t)
        
        Tx=reshape(t(1:nxx),[endo_nbr,endo_nbr,regime_nbr]);
        
        Tsig=reshape(t(nxx+(1:endo_nbr*regime_nbr)),[endo_nbr,1,regime_nbr]);
        
        Te=reshape(t(nxx+endo_nbr*regime_nbr+1:end),[endo_nbr,exo_nbr,regime_nbr]);
        
        ss1=ss;
        
    end

    function [t1,f]=iterator(t0)
        
        TT=struct();
        
        [TT.Tx,TT.Tsig,TT.Te,TT.ss]=resizer(t0);
        
        [mat,Resids,Q,c,own_stuff]=subutils.load_state_matrices(obj,structural_matrices,TT,x0,inner_probs,control_shocks,own_stuff);
        
        T1 = TT;
        
        for ii=1:regime_nbr
            
            pijAjTj=0;
            
            for jj=1:regime_nbr
                
                pijAjTj = pijAjTj + mat.Aplus(:,:,ii,jj) * TT.Tx(:,:,jj);
                
            end
            
            tmp=(pijAjTj + mat.A0(:,:,ii));
            
            if ii==1
                
                U=tmp(:,:,ones(1,regime_nbr));
                
            else
                
                U(:,:,ii)=tmp;
                
            end
            
            tmp=-tmp \ aie;
            
            T1.Tx(:,:,ii) = tmp * mat.Aminus(:,:,ii);
            
            T1.Te(:,:,ii) = tmp * mat.B(:,:,ii);
            
        end
        
        % non-certainty equivalence
        %--------------------------
        if any(Resids(:))
            
            biga=zeros(endo_nbr*regime_nbr);
            
            for s0=1:regime_nbr
                
                rrr=endo_nbr*(s0-1)+(1:endo_nbr);
                
                for s1=1:regime_nbr
                    
                    ccc=endo_nbr*(s1-1)+(1:endo_nbr);
                    
                    biga(rrr,ccc)=mat.Aplus(:,:,s0,s1);
                    
                    if s0==s1
                        
                        biga(rrr,ccc)=biga(rrr,ccc)+U(:,:,s0);
                        
                    end
                    
                end
                
            end
            
            T1.Tsig=reshape(biga\Resids(:),[endo_nbr,1,regime_nbr]);
            
        else
            
            T1.Tsig=zeros(endo_nbr,1,regime_nbr);
            
        end
        
        t1 = subutils.vectorizer(T1);
        
        f = abs(t1 - t0);
        
    end

end