function [welf,retcode,V,d]=calculate_loss(m,lossstr,shocksdb)
% Calculates welfare
%
% ::
%
%   [welf,retcode,V,d]=calculate_loss(m,lossstr,shocksdb)
%
% Args:
%
%    - **m** (struct | rise | dsge): scalar or vector of model objects.
%
%    - **lossstr** (char): loss function
%
%    - **shocksdb** (ts | struct | empty): Shocks database to condition on over
%      simulation. if empty, the unconditional welfare is returned. if not
%      empty the conditional welfare is calculated instead
%
%    - **varargin** : optional arguments coming in pairs
%
% Returns:
%
%    - **welf** [scalar|vector]: conditional loss (scalar) or unconditional
%      loss (scalar if number of regimes = 1)
%
%    - **retcode** [scalar]: return code
%
%    - **V** [n x n x h]: array of equilibrium matrices
%
%    - **d** [h x 1 vector]: constant in loss
%
% Note:
%
%    The loss is such that :math:`v(x_t,r_t)=x_t.'*Vrt*x_t+drt`
%

callers=dbstack;

callers={callers.name};

if any(strcmp(callers,'estimate'))
    
    retcode=0;
    
else
    
    [m,retcode]=dsge_tools.solve.solve(m);
    
end

welf=[]; V=[]; d=[];

if all(retcode)

    return

end

fslv=~true; % just for debugging purposes one can use fsolve instead of tfqmr

npar=dsge_tools.get.number_of_parameterizations(m);
    
V=cell(1,npar);

welf=V;

d=V;

for ip=1:npar

    [welf{ip},V{ip},d{ip}]=loss_engine(...
        dsge_tools.solve.load_state_space(m,ip),lossstr,shocksdb,fslv,retcode(ip));

end

welf=cell2mat(welf);

V=cell2mat(V);

d=cell2mat(d);

end


function [welf,V,d]=loss_engine(m,lossstr,shocksdb,fslv,retcode)

welf=[]; V=[]; d=[];

if retcode

    return

end

[Tx,Te]=dsge.set_solution_to_companion(m);

[n,n]=size(Tx{1});

h=numel(Tx);

ss=dsge_tools.get.stst_(m);

lv=dsge_tools.get.log_vars(m);

ss(lv,:)=exp(ss(lv,:));

h=numel(Tx);

weightCell=dsge.compile_loss_weights(m,lossstr);

beta=weightCell{1}.discount;

param=dsge_tools.get.param_(m);

npcols=size(param,2);

W0=full(load_weight(ss(:,1),param(:,1)));

for ih=2:h
    
    jh=min(ih,npcols);
    
    W0(:,:,ih)=load_weight(ss(:,jh),param(:,jh));
    
end

Q=m.state_space{1}.transition_matrices.Q;

if size(W0,3)==1 && h>1
    
    W0=W0(:,:,ones(1,h));
    
end

if numel(beta)>1
    
    if ~all(beta==beta(1))
        
        error(['The formulae need to be modified to handle the case of ',...
            'different discounts. Send an email to junior.maih@gmail.com'])
        
    end
    
    beta=beta(1);
    
end

ping=find(any(W0(:,:,1)));

for ih=2:h

    ping=union(ping,find(any(W0(:,:,ih))));

end

% PATCH 1 (value-recursion truncation). The recursion below propagates V
% through Tx(ping,ping). Restricting to the loss-relevant variables is exact
% only if they form a closed subsystem under Tx, which they generally do not
% (they depend on the shocks, the instrument and, under optimal policy, the
% Lagrange multipliers). Run the recursion on the FULL state vector instead;
% W0 already has zeros outside the loss variables, so this is the "whole T,
% W padded with zeros" formulation. Non-loss variables that are also non-state
% have zero rows in W and zero columns in Tx and contribute nothing, so the
% result is exact -- verified to match an independent moment calculation.
ping=1:n;

np=numel(ping);

W=W0(ping,ping,:);

V0 = W; % V0=zeros(np,np,h);

if beta>0
    
    if fslv
        
        options=struct('Display','iter',...
            'TolFun',m.options.fix_point_TolFun,...
            'MaxIter',m.options.fix_point_maxiter);
        
        [V,fval,exitflag]=fsolve(@fslvObjective,V0,options);
        
        retcode=utils.optim.exitflag(exitflag,V,fval,m.options.fix_point_TolFun);
        
    else
        
        B=W(:);
        
        [V,retcode]=transpose_free_quasi_minimum_residual(@objective,... % coefficient matrix
            B(:),... % right hand side
            V0(:),... % initial guess
            m.options.fix_point_TolFun,... % tolerance level
            m.options.fix_point_maxiter,... % maximum number of iterations
            m.options.fix_point_verbose); % flag for printing progress or not

        % tfqmr can stop short of its tolerance (it exits on "stagnation"
        % with retcode 0). A partially solved V makes the loss jagged in
        % the policy parameters and hands optimizers spurious low pockets,
        % so verify the residual and solve exactly if needed.
        [V,retcode]=verify_or_solve_exactly(V,retcode,B,np,h,Q,Tx,ping,...
            beta,m.options.fix_point_TolFun,@objective);

    end
    
else
    
    retcode=0;
    
    V=V0;
    
end

if retcode
    
    return % error(decipher(retcode))
    
end

V=reshape(full(V),[np,np,h]);

d=zeros(h,1);

if beta>0
    
    G=zeros(h,1);
    
    for iii=1:h
        
        G(iii)=trace(V(:,:,iii)*(Te{iii}(ping,:)*Te{iii}(ping,:).'));
        
    end
    
    % PATCH 2 (operator precedence). '*' and '\' associate left to right, so
    % the unbracketed  beta*(eye(h)-beta*Q)\Q*G  is ((beta*(eye(h)-beta*Q))\Q)*G
    % = (1/beta)*(I-beta*Q)^{-1}*Q*G -- the constant is too large by 1/beta^2.
    d=beta*((eye(h)-beta*Q)\(Q*G));
    
end

simulator=dsge.compile_conditioning_shocks(m,shocksdb);


[simuls,retcode]=simulator(m);

welf = -10e+8;

if retcode
    
    return
    
    
end

if isempty(simuls)
    
    [welf1,retcode]=unconditional_welfare();
    
else
    
    [welf1,retcode]=conditional_welfare();
    
end

if ~retcode
    
    welf = welf1;
    
end


    function [welf,retcode]=conditional_welfare()
        
        retcode=0;
        
        simuls=pages2struct(simuls);
        
        endo_names=m.endogenous.name;
        
        span=get(simuls.(endo_names{1}),'NumberOfObservations');
        
        x=nan(n,span);
        
        for iv=1:n
            
            x(iv,:)=double(simuls.(endo_names{iv}));
            
        end

        if isfield(simuls,'regime')

            regimes=double(simuls.regime);

            Vt=V;

            dt=d;

            sst=ss;

        else

            regimes=ones(1,span);

            Vt=mean(V,3);

            dt=mean(d);

            sst=mean(ss,2);

        end

        
        welf=0;        
        
        for t=1:span
            
            rt=regimes(t);
            
            if isnan(rt)
                
                continue
                
            end
            
            xt=x(ping,t)-sst(ping,rt);
            
            welf=welf+xt.'*Vt(:,:,rt)*xt+dt(rt);
            
        end
        
        welf=welf/span;

    end


    function [welf,retcode]=unconditional_welfare()

        m.options.autocov_model_resolve=false;
        
        [ACV,retcode]=dsge_tools.moments.theoretical_autocovariances(m,0,ping);
        
        if retcode
            
            welf=[];
            
            return
            
        end
        
        welf=nan(1,h);
        
        for ii=1:h
            
            welf(ii)=d(ii)+trace(V(:,:,ii)*ACV{ii});
            
            % dynare style
            % wi=W(:,:,ii); acvi=ACV(:,:,ii); welf(ii) = full(wi(:).'*acvi(:));
            
        end
        
        %    EVV=sum(diag(W0{1}).*diag(ACV))/(1-beta{1});
        
        % EV == EVV?
    end


    function ax=fslvObjective(V)
        
        ax=objective(V);
        
        ax=reshape(ax,[np,np,h])-W;
        
    end


    function ax=objective(V)
        
        V=reshape(full(V),[np,np,h]);
        
        F=zeros(np,np,h);
        
        for ii=1:h
            
            for jj=1:h
                
                F(:,:,ii)=F(:,:,ii)+...
                    Q(ii,jj)*Tx{ii}(ping,ping).'*V(:,:,jj)*Tx{ii}(ping,ping);
                
            end
            
            F(:,:,ii)=V(:,:,ii)-beta*F(:,:,ii);
            
        end
        
        ax=F(:);
        
    end


    function w=load_weight(ss,param)
                        
        w=utils.code.evaluate_functions(weightCell,ss,param);
        
        w=reshape(w,n,n);

    end

end


function [V,retcode]=verify_or_solve_exactly(V,retcode,B,np,h,Q,Tx,ping,beta,tol,objective)
% Accepts the iterative solution of
%
%   V_i - beta*sum_j Q(i,j)*T_i'*V_j*T_i = W_i,   T_i=Tx{i}(ping,ping)
%
% only if it converged and its relative residual is within tol. Otherwise
% solves the vec form of the same system exactly,
%
%   (I - beta*[Q(i,j)*kron(T_i',T_i')]_{ij}) vec(V) = vec(W),
%
% with a sparse direct solver when the matrix is affordable, and returns a
% nonzero retcode (201: no convergence) when it is not (or when the direct
% solve fails).

max_nnz=1e7; % ~160MB of sparse storage for the coefficient matrix

if isempty(tol)

    tol=sqrt(eps);

end

tol=max(tol,10*eps);

nB=max(norm(B),realmin);

if ~retcode && all(isfinite(V))

    relres=norm(objective(V)-B)/nB;

    if relres<=tol

        return

    end

end

N=np^2*h;

T=cell(1,h);

nnzA=N;

for ii=1:h

    T{ii}=sparse(Tx{ii}(ping,ping)).';

    nnzA=nnzA+nnz(T{ii})^2*nnz(Q(ii,:));

end

failed=201;

if nnzA>max_nnz

    retcode=failed;

    return

end

blocks=cell(h);

for ii=1:h

    K=kron(T{ii},T{ii});

    for jj=1:h

        blocks{ii,jj}=(-beta*Q(ii,jj))*K;

    end

end

A=speye(N)+cell2mat(blocks);

ws=warning('off','MATLAB:singularMatrix');

ws2=warning('off','MATLAB:nearlySingularMatrix');

Vx=A\B;

warning(ws2); warning(ws);

if all(isfinite(Vx)) && norm(A*Vx-B)/nB<=max(tol,sqrt(eps))

    V=Vx;

    retcode=0;

else

    retcode=failed;

end

end
