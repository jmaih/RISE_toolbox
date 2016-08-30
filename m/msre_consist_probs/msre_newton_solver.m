function [T,retcode,c,Q]=msre_newton_solver(obj,inner_probs,T0,x0,control_shocks,options,kron_method)

if nargin < 7
    
    kron_method = true;
    
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

end

n=obj.endogenous.number;

exo_nbr=sum(obj.exogenous.number);

h=obj.markov_chains.regimes_number;

if isempty(T0);
    
    T0.Tx=zeros(n,n,h);
    
    T0.Tsig=zeros(n,1,h);
    
    T0.Te=zeros(n,exo_nbr,h);
    
    T0.ss=[];
    
end

if isempty(control_shocks)
    
    control_shocks=randn(exo_nbr,1,10000);
    
end

order_var=obj.order_var;

inner_probs=subutils.process_probabilities(obj,inner_probs);

aie=speye(n);

% Place holders for constant state matrices
own_stuff=[];

nxx=numel(T0.Tx);

if ~isfield(options,'local_approximation')
    
    options.local_approximation=false;
    
end

if options.local_approximation
    
    obj.routines.steady_state_model=vector_to_steady_state_model(x0);
    
    obj=set(obj,'steady_state_imposed',true);
    
end

[obj,structural_matrices,retcode]=compute_steady_state(obj);

if retcode
    
    error(decipher(retcode))
    
end

ss=obj.solution.ss;
    
for ireg=1:numel(ss)
    
    ss{ireg}=ss{ireg}(order_var);
    
end

x0=x0(order_var);

n_npb=n*n;

if kron_method
    
    G=zeros(n_npb*h);
    
else
    
    LMINUS=cell(1,h);
    
    LPLUS=cell(h);
    
end

I_nx_nd=speye(n_npb);

best_solution=struct('f',inf);

[t,~,retcode]=fix_point_iterator(@iterator,subutils.vectorizer(T0),options);

if ~isfield(options,'use_best') 
    
    options.use_best=false;
    
end

if options.use_best
    
    if retcode
        
        warning('Using best found')
        
        retcode=0;
        
    end
    
    T=best_solution.T;
    
    c=best_solution.c;
    
    Q=best_solution.Q;
    
else
    
    T=struct();
    
    [T.Tx,T.Tsig,T.Te,T.ss]=resizer(t);
    
end

    function [Tx,Tsig,Te,ss1]=resizer(t)
        
        Tx=reshape(t(1:nxx),[n,n,h]);
        
        Tsig=reshape(t(nxx+(1:n*h)),[n,1,h]);
        
        Te=reshape(t(nxx+n*h+1:end),[n,exo_nbr,h]);
        
        ss1=ss;
        
    end

    function [t1,f]=iterator(t0)
        
        TT=struct();
        
        [TT.Tx,TT.Tsig,TT.Te,TT.ss]=resizer(t0);
        
        [mat,Resids,Q,c,own_stuff]=subutils.load_state_matrices(obj,structural_matrices,TT,x0,inner_probs,control_shocks,own_stuff);
        
        T1 = TT;
        
        W=TT.Tx;
        
        for r0=1:h
            
            U=mat.A0(:,:,r0);
            
            for r1=1:h
                
                U = U + mat.Aplus(:,:,r0,r1) * TT.Tx(:,:,r1);
                
            end
            
            Ui=U\aie;
            
            if r0==1
                
                UUUUUU=U(:,:,ones(1,h));
                
            else
                
                UUUUUU(:,:,r0)=U;
                
            end
            
            T1_fi=-Ui*mat.Aminus(:,:,r0);
            
            W(:,:,r0)=W(:,:,r0)-T1_fi;
            
            Lminus=-T1_fi;
            
            if kron_method
                
                Lminus_prime=Lminus.';
                
            else
                
                LMINUS{r0}=sparse(Lminus);
                
            end
            
            rows=(r0-1)*n_npb+1:r0*n_npb;
            
            for r1=1:h
                
                cols=(r1-1)*n_npb+1:r1*n_npb;
                
                Lplus01=Ui*mat.Aplus(:,:,r0,r1);
                
                if kron_method
                    % build G
                    %--------
                    tmp=kron(Lminus_prime,Lplus01);
                    
                    if r0==r1
                        
                        tmp=tmp-I_nx_nd;
                        
                    end
                    
                    G(rows,cols)=tmp;
                    
                else
                    
                    LPLUS{r0,r1}=sparse(Lplus01);
                    
                end
                
            end
            
            
        end
        
        W=reshape(W,[n,n*h]);
        
        if kron_method
            % update T
            %---------
            delta=G\W(:);
            
        else
            
            delta0=[];
            
            [delta,retcode]=utils.optim.linear_systems_solver(...
                @(x)find_newton_step(x,LPLUS,LMINUS),-W(:),delta0,options);
        end
        
        T1.Tx=TT.Tx+reshape(delta,[n,n,h]);
        
        update_the_rest()
            
        % non-certainty equivalence
        %--------------------------
        
        t1 = subutils.vectorizer(T1);
        
        f = abs(t1 - t0);
        
        max_f = max(f);
        
        if max_f < best_solution.f
            
            best_solution.f = max_f;
            
            best_solution.T=T1;
            
            best_solution.c=c;
            
            best_solution.Q=Q;
            
        end
        
        function update_the_rest()
            
%             W1=T1.Tx;
            
            is_certainty_equivalent = ~any(Resids(:));
            
            if is_certainty_equivalent
                
                T1.Tsig=zeros(n,1,h);
                
            else
                
                biga=zeros(n*h);
                
            end
            
            for r00=1:h
                
                U=mat.A0(:,:,r00);
                
                for r11=1:h
                    
                    U=U+mat.Aplus(:,:,r00,r11)*T1.Tx(:,:,r11);
                    
                end
                
                Ui=U\aie;
                
                if ~is_certainty_equivalent
                    
                    rrr=n*(r00-1)+(1:n);
                    
                    for r11=1:h
                        
                        ccc=n*(r11-1)+(1:n);
                        
                        biga(rrr,ccc)=mat.Aplus(:,:,r00,r11);
                        
                        if r00==r11
                            
                            biga(rrr,ccc)=biga(rrr,ccc)+U;
                            
                        end
                        
                    end
                    
                end
                
%                 new_T1_fi = -Ui*mat.Aminus(:,:,r00);
                
%                 W1(:,:,r00)=W1(:,:,r00)-new_T1_fi;
                
                T1.Te(:,:,r00) = -Ui * mat.B(:,:,r00);
                
            end
            
            if ~is_certainty_equivalent
                
                T1.Tsig=reshape(biga\Resids(:),[n,1,h]);
                
            end
            
        end
        
        function Gd=find_newton_step(delta,Lplus,Lminus)
            
            Gd=zeros(n*n,h); % G*delta
            
            delta=reshape(delta,[n*n,h]);
            
            for r00=1:h
                
                for r11=1:h
                    
                    Gd(:,r00)=Gd(:,r00)+vec(Lplus{r00,r11}*reshape(delta(:,r11),n,n)*Lminus{r00});
                    
                end
                
            end
            
            Gd=delta-Gd;
            
            Gd=Gd(:);
            
        end
        
    end

end

function ssm=vector_to_steady_state_model(y0)

ssm=@numerical_ssmodel;

    function [yout,param]=numerical_ssmodel(~,~,~,param,~,~,~)
        
        yout=y0;
        
    end

end