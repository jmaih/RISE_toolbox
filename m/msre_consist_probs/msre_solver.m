function [T,retcode,c,Q]=msre_solver(obj,inner_probs,T0,x0,control_shocks,options)

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

process_probabilities()

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
        
[t,~,retcode]=fix_point_iterator(@iterator,vectorizer(T0),...
    options);

T=struct();

[T.Tx,T.Tsig,T.Te,T.ss]=resizer(t);

    function process_probabilities()
        
        regimes=cell2mat(obj.markov_chains.regimes(2:end,2:end));
        
        ordered_names=obj.endogenous.name(obj.order_var);
        
        do_replace=@replacement; %#ok<NASGU>
        
        express=['\<',parser.cell2matize(obj.endogenous.name),'\>'];
        
        for iprob=1:size(inner_probs,1)
            
            prob=regexprep(inner_probs{iprob,2},express,'${do_replace($1)}');
            
            [~,markovChainLocs]=decompose_name(inner_probs{iprob,1});
            
            disp(prob)
            
            inner_probs{iprob,2}=struct('func',str2func(['@(x)',prob]),...
                'chain_loc',markovChainLocs);
            
        end
        
        function out=replacement(vname)
            
            ploc=int2str(find(strcmp(vname,ordered_names)));
            
            out=['x(',ploc,',:)'];
            
        end
        
        function [currState,markovChainLocs]=decompose_name(tp_name)
            
            underscore=find(tp_name=='_');
            
            currState=str2double(tp_name(underscore(2)+1:underscore(3)-1));
            
            chain_name=tp_name(1:underscore(1)-1);
            
            chain_loc= strcmp(chain_name,obj.markov_chains.chain_names);
            
            markovChainLocs = regimes(:,chain_loc) == currState;
            
        end
        
    end

    function [Tx,Tsig,Te,ss1]=resizer(t)
        
        Tx=reshape(t(1:nxx),[endo_nbr,endo_nbr,regime_nbr]);
        
        Tsig=reshape(t(nxx+(1:endo_nbr*regime_nbr)),[endo_nbr,1,regime_nbr]);
        
        Te=reshape(t(nxx+endo_nbr*regime_nbr+1:end),[endo_nbr,exo_nbr,regime_nbr]);
        
        ss1=ss;
        
    end

    function t=vectorizer(T)
        
        t=[T.Tx(:);T.Tsig(:);T.Te(:)];
        
    end

    function [t1,f]=iterator(t0)
        
        TT=struct();
        
        [TT.Tx,TT.Tsig,TT.Te,TT.ss]=resizer(t0);
        
        [mat,Resids,Q,c,own_stuff]=load_state_matrices(obj,structural_matrices,TT,x0,inner_probs,control_shocks,own_stuff);
        
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
                
        t1 = vectorizer(T1);
        
        f = abs(t1 - t0);
        
    end

end

function [mat,Resids,Q,c,own_stuff]=load_state_matrices(m,structural_matrices,T,x0,inner_probabilities,control_shocks,own_stuff)

if nargin<6
    % this is to avoid using persistent variables
    own_stuff=[];
    
end

if isempty(own_stuff)
    
    own_stuff=struct();
    
    own_stuff.is_constant = isempty(inner_probabilities);
    
    own_stuff.derivs=struct();
    
end

if own_stuff.is_constant
    
    c=[];
    
    if ~isempty(own_stuff.derivs.A00)
        
        mat.Aplus=own_stuff.derivs.Aplus0;
        
        mat.A0=own_stuff.derivs.A00;
        
        mat.Aminus=own_stuff.derivs.Aminus0;
        
        mat.B=own_stuff.derivs.B0;
        
        Resids=own_stuff.Resids;
        
        Q=own_stuff.Q;
        
        return
        
    end
    
else
    
    c=simulated_probabilities(T,x0,inner_probabilities,control_shocks);
    
    m=set(m,'parameters',[inner_probabilities(:,1),num2cell(c)']);
    
end

[structural_matrices,retcode]=dsge_tools.evaluate_all_derivatives(m,structural_matrices,[]);

[pos,siz,shock_horizon]=dsge_tools.rehash_topology(m,structural_matrices);

[sm]=utils.solve.pull_first_order_partitions(structural_matrices.dv,pos.v);

Resids=structural_matrices.user_resids;

adjusted=struct();
adjusted.bf_cols=pos.t.bf;
adjusted.pb_cols=pos.t.pb;

adjusted.nd=siz.nd;
adjusted.siz=siz; % adjusted sizes
accelerate=false;
% accelerate=options.solve_accelerate && siz.ns;

[sm,adjusted,accelSupport]=dsge_tools.aggregate_matrices(sm,siz,adjusted,accelerate);

[mat.Aplus,mat.A0,mat.Aminus]=dsge_tools.full_state_matrices(siz,sm);

Q = structural_matrices.transition_matrices.Q;

mat.B=full(sm.de_0{1});

mat.B=mat.B(:,:,ones(siz.h,1));

for istate = 2:siz.h
    
    mat.B(:,:,istate)=sm.de_0{istate};
    
end

if own_stuff.is_constant
    
    own_stuff.derivs.Aplus0=mat.Aplus;
    
    own_stuff.derivs.A00=mat.A0;
    
    own_stuff.derivs.Aminus0=mat.Aminus;
    
    own_stuff.derivs.B0=mat.B;
    
    own_stuff.Resids=Resids;
    
    own_stuff.Q=Q;
    
end

end

