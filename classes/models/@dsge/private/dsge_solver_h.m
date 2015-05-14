function [T,eigval,retcode,obj]=dsge_solver_h(obj,structural_matrices)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

% the obj going out probably contains the changed options

if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    T=struct();
    return
end

% It is assumed the steady states have been solved and the derivatives
% evaluated at different orders
%%
Pfunc=@utils.kronecker.sum;

%% begin
T=struct();
% options.solve_order 1
%--------
if obj.options.solve_order>=1
    pos=obj.locations.before_solve;
    siz=obj.siz.before_solve;
    % collect the sizes
    %------------------
    siz.h=size(structural_matrices.dv,2);
    siz.nT=siz.ns+siz.np+siz.nb+siz.nf;
    shock_horizon=max(obj.exogenous.shock_horizon);
    siz.nz=siz.np+siz.nb+1+siz.ne*(1+shock_horizon);
    siz.nd=size(structural_matrices.dv{1,1},1); % number of equations
    pos.z.e_plus=pos.z.e_0(end)+(1:shock_horizon*siz.ne);
    
    % Structure of elements that will move across different orders
    %-------------------------------------------------------------
    others=struct();
    others.theta_hat=structural_matrices.theta_hat;
    
    [T.Tz,others,eigval,retcode,obj.options]=solve_first_order(structural_matrices,...
        others,siz,pos,obj.options,shock_horizon);
    
    % higher orders
    %--------------
    if obj.options.solve_order>=2 && ~retcode
        [T,retcode]=solve_higher_orders(T,others,obj.options.solve_accelerate);
    end
end

    function [T,retcode]=solve_higher_orders(T,others,accelerate)
        % second-order moments
        %---------------------
        [mompos,zz_sig_loc]=moments_positions(siz,2);
        Eu_u=sparse(siz.nz^2,siz.nz^2);
        Eu_u(mompos,zz_sig_loc)=1;
        
        a0_z=zeros(siz.nv,siz.nz);
        a0_z(pos.v.b_minus,pos.z.b)=eye(siz.nb);
        a0_z(pos.v.p_minus,pos.z.p)=eye(siz.np);
        a0_z(pos.v.e_0,pos.z.e_0)=eye(siz.ne);
        a1_z=zeros(siz.nv,siz.nz);
        
        hz=zeros(siz.nz);
        hz(pos.z.sig,pos.z.sig)=1;
        hz(siz.np+siz.nb+1+(1:siz.ne*shock_horizon),pos.z.e_plus)=eye(shock_horizon*siz.ne);
        dbf_plus=others.dbf_plus;
        for rt=1:siz.h
            for rplus=1:siz.h
                % this preconditioning could also be done once and for all and
                % should be included in Aplus
                dbf_plus{rt,rplus}=others.Ui(:,:,rt)*dbf_plus{rt,rplus};
            end
        end
        
        Dzz=second_order_rhs();
        [T.Tzz,retcode]=solve_generalized_sylvester(Dzz,2);
        clear Dzz
        
        if obj.options.solve_order>2 && ~retcode
            % store a0, a1 and hz as sparse since we know their non-zero
            % pattern (most likely) will not change.
            a0_z=sparse(a0_z);
            a1_z=sparse(a1_z);
            hz=sparse(hz);
            % third-order moments: no skewed shocks
            %--------------------------------------
            % Eu_u_u=sparse(siz.nz^3,siz.nz^3);
            
            Dzzz=third_order_rhs();
            [T.Tzzz,retcode]=solve_generalized_sylvester(Dzzz,3);
            clear Dzzz
            
            if obj.options.solve_order>3 && ~retcode
                % fourth-order moments
                %---------------------
                % Eu_u_u_u=sparse(siz.nz^4,siz.nz^4);
                Dzzzz=fourth_order_rhs();
                [T.Tzzzz,retcode]=solve_generalized_sylvester(Dzzzz,4);
                clear Dzzzz
                
                if obj.options.solve_order>4 && ~retcode
                    Dzzzzz=fifth_order_rhs();
                    [T.Tzzzzz,retcode]=solve_generalized_sylvester(Dzzzzz,3);
                    clear Dzzzzz
                    
                    if obj.options.solve_order>5
                        error('perturbations of order greater than 5 not implemented');
                    end
                end
            end
        end
        
        function Dzz=second_order_rhs()
            Dzz=preallocate_rhs(2);
            for r0=1:siz.h
                a0_z(pos.v.t_0,:)=T.Tz{r0};
                hz(pos.z.pb,:)=T.Tz{r0}(pos.t.pb,:);
                for r1=1:siz.h
                    a1_z(pos.v.bf_plus,:)=T.Tz{r1}(pos.t.bf,:);
                    a0_z(pos.v.bf_plus,:)=T.Tz{r1}(pos.t.bf,:)*hz;
                    a0_z(pos.v.theta_plus,pos.z.sig)=others.theta_hat{r0,r1};
                    Dzz(:,:,r0)=Dzz(:,:,r0)+dvv_Evz_vz();
                end
                % precondition
                Dzz(:,:,r0)=-others.Ui(:,:,r0)*Dzz(:,:,r0);
            end
            
            function res=dvv_Evz_vz()
                res=fvv_vx_vx(structural_matrices.dvv{r0,r1},a0_z);
                res=res+fvv_vx_vx(structural_matrices.dvv{r0,r1},a1_z)*Eu_u;
                if obj.options.debug
                    tmp1=fvv_vx_vx(structural_matrices.dvv{r0,r1},a1_z)*Eu_u;
                    tmp2=structural_matrices.dvv{r0,r1}*kron(a1_z,a1_z)*Eu_u;
                    max(max(abs(tmp1-tmp2)))
                    keyboard
                end
            end
        end
        
        function Dzzz=third_order_rhs()
            Dzzz=preallocate_rhs(3);
            a0_zz=zeros(siz.nv,siz.nz^2);
            a1_zz=zeros(siz.nv,siz.nz^2);
            hzz=zeros(siz.nz,siz.nz^2);
            for r0=1:siz.h
                a0_z(pos.v.t_0,:)=T.Tz{r0};
                a0_zz(pos.v.t_0,:)=T.Tzz{r0};
                hz(pos.z.pb,:)=T.Tz{r0}(pos.t.pb,:);
                hzz(pos.z.pb,:)=T.Tzz{r0}(pos.t.pb,:);
                for r1=1:siz.h
                    a1_z(pos.v.bf_plus,:)=T.Tz{r1}(pos.t.bf,:);
                    a0_z(pos.v.bf_plus,:)=T.Tz{r1}(pos.t.bf,:)*hz;
                    a0_z(pos.v.theta_plus,pos.z.sig)=others.theta_hat{r0,r1};
                    
                    a0_zz(pos.v.bf_plus,:)=T.Tz{r1}(pos.t.bf,:)*hzz+T.Tzz{r1}(pos.t.bf,:)*kron(hz,hz);
                    a1_zz(pos.v.bf_plus,:)=T.Tzz{r1}(pos.t.bf,:);
                    
                    Dzzz(:,:,r0)=Dzzz(:,:,r0)+...
                        dvvv_Evz_vz_vz()+...
                        dvv_Evz_vzz()+...
                        others.dbf_plus{rt,rplus}*Tzz_hz_hzz();
                end
                % precondition
                Dzzz(:,:,r0)=-others.Ui(:,:,r0)*Dzzz(:,:,r0);
            end
            
            function res=dvv_Evz_vzz()
                Eu_u_hz=kron(Eu_u,hz);
                Eu_hz_u=utils.kronecker.shuffle(Eu_u_hz,[siz.nz,siz.nz],[siz.nz,siz.nz].^2);
                res=fvv_vx_vxx_omega_1(structural_matrices.dvv{r0,r1},a0_z,a0_zz+a1_zz*Eu_u)+...
                    structural_matrices.dvv{r0,r1}*kron(a1_z,a1_zz)*(Eu_u_hz+Eu_hz_u);
            end
            
            function res=Tzz_hz_hzz()
                res=fvv_vx_vxx_omega_1(T.Tzz{r1}(pos.t.bf,:),hz,hzz);
            end
            
            function res=dvvv_Evz_vz_vz()
                res=fvvv_vx_vx_vx(structural_matrices.dvvv{r0,r1},a0_z);
                k011u=kron(a0_z,kron(a1_z,a1_z)*Eu_u);
                siz_a0=size(a0_z);
                siz_a1_u=siz_a0;
                k11u0 = utils.kronecker.shuffle(k011u,siz_a0,siz_a1_u.^2);
                k1u01 = utils.kronecker.shuffle(k011u,siz_a0.*siz_a1_u,siz_a1_u);
                res=res+structural_matrices.dvvv{r0,r1}*(k011u+k11u0+k1u01);
            end
        end
        
        function D=preallocate_rhs(oo)
            D=zeros(siz.nd,siz.nz^oo,siz.h);
        end
        
        function [X,retcode]=solve_generalized_sylvester(D,oo)
            % [X,retcode] = tfqmr(@afun,D(:),obj.options.fix_point_TolFun);
            %     0 tfqmr converged to the desired tolerance TOL within MAXIT iterations.
            %     1 tfqmr iterated MAXIT times but did not converge.
            %     2 preconditioner M was ill-conditioned.
            %     3 tfqmr stagnated (two consecutive iterates were the same).
            %     4 one of the scalar quantities calculated during tfqmr became too
            
            % shrink D and so on
            %--------------------
            if accelerate
                [shrink,expand]=utils.kronecker.shrink_expand(siz.nz,oo);
                nkept=sum(shrink);
                if obj.options.debug
                    Dtest=D;
                end
                D=D(:,shrink,:);
                if obj.options.debug
                    Dbig=D(:,expand,:);
                    max(abs(Dtest(:)-Dbig(:)))
                end
            else
                if obj.options.debug
                    [shrink,expand]=utils.kronecker.shrink_expand(siz.nz,oo);
                    D=D(:,shrink,:);
                    D=D(:,expand,:);
                end
                nkept=siz.nz^oo;
            end
            % compute E(kron(u,...,u))
            E_u_o=Ekron_u_oo();
            
            % expand final result: maybe, maybe not
            %--------------------------------------
            [X,retcode]=transpose_free_quasi_minimum_residual(@afun,D(:),... % right hand side
                [],... %x0 initial guess
                obj.options.fix_point_TolFun,... % tolerance level
                obj.options.fix_point_maxiter,... % maximum number of iterations
                obj.options.fix_point_verbose);
            if ~retcode
                if accelerate
                    X=reshape(X,[siz.nd,nkept,siz.h]);
                    X=X(:,expand,:);
                else
                    X=reshape(X,[siz.nd,siz.nz^oo,siz.h]);
                end
                tmp=cell(1,siz.h);
                for r0=1:siz.h
                    tmp{r0}=X(:,:,r0);
                end
                X=tmp;
            end
            
            function ATC_plus_T=afun(tau)
                tau=reshape(tau,[siz.nd,nkept,siz.h]);
                ATC_plus_T=zeros(size(tau));
                for r00=1:siz.h
                    hz(pos.z.pb,:)=T.Tz{r00}(pos.t.pb,:);
                    AT=0;
                    for r1=1:siz.h
                        AT=AT+dbf_plus{r00,r1}*tau(pos.t.bf,:,r1);
                    end
                    if accelerate
                        AT=AT(:,expand);
                    end
                    ATCzzz=temporary_term();
                    ATC_plus_T(:,:,r00)=ATCzzz+tau(:,:,r00);
                end
                ATC_plus_T=ATC_plus_T(:);
                function ATCzzz=temporary_term()
                    % core terms
                    %-----------
                    ATCzzz=utils.kronecker.A_times_k_kron_B(AT,hz,oo);
                    ATCzzz=ATCzzz+AT*E_u_o;
                    if oo==3
                        ATCzzz=ATCzzz+AT*Pfunc(hz,Eu_u);
                    end
                    if accelerate
                        ATCzzz=ATCzzz(:,shrink);
                    end
                end
            end
            function E_u_o=Ekron_u_oo()
                switch oo
                    case 2
                        E_u_o=Eu_u;
                    case 3
                        % no skewed shocks?
                        E_u_o=0;
                    otherwise
                        error(['approximation of order ',int2str(oo),' not yet implemented'])
                end
            end
        end
    end
end

function Y=fvv_vx_vx(dvv,vz,tol)
if nargin<3
    tol=sqrt(eps);
end
if max(abs(dvv(:)))<tol || max(abs(vz(:)))<tol
    nd=size(dvv,1);
    [~,nz]=size(vz);
    Y=zeros(nd,nz^2);
else
    Y=utils.kronecker.A_times_k_kron_B(dvv,vz,2);
end

end

function Y=fvvv_vx_vx_vx(dvvv,vz,tol,new_algo)
if nargin<4
    new_algo=true;
    if nargin<3
        tol=sqrt(eps);
    end
end
[nd]=size(dvvv,1);
[nv,nz]=size(vz);

if max(abs(dvvv(:)))<tol || max(abs(vz(:)))<tol
    Y=zeros(nd,nz^3);
else
    if new_algo
        Y=utils.kronecker.A_times_k_kron_B(dvvv,vz,3);
    else
        Y=reshape((reshape(dvvv,nd*nv^2,nv)*vz)',nz*nd*nv,nv);
        Y=reshape((Y*vz)',nz*nz*nd,nv);
        Y=permute(reshape(full(Y*vz),nz,nz,nd,nz),[3,2,1,4]);
        Y=reshape(Y,[nd,nz^3]);
    end
end

end

function [Y]=fvv_vx_vxx_omega_1(dvv,vz,vzz,tol,new_algo)
if nargin<5
    new_algo=true;
    if nargin<4
        tol=sqrt(eps);
    end
end
[nd]=size(dvv,1);
[~,nz]=size(vz);
kk=2;ll=3;jj=4;
if max(abs(dvv(:)))<tol || max(abs(vz(:)))<tol || max(abs(vzz(:)))<tol
    Y=zeros(nd,nz^3);
else
    if new_algo
        Y=utils.kronecker.A_times_B_kron_C(dvv,vz,vzz);
        Y=reshape(Y,[nd,nz,nz,nz]);
        Y=Y+ipermute(Y,[1,jj,ll,kk])+ipermute(Y,[1,jj,kk,ll]);
    else
        Y=reshape((reshape(dvv,nd*nv,nv)*vzz)',nz^2*nd,nv);
        Y=permute(reshape(full(Y*vz),nz,nz,nd,nz),[ll,1,kk,jj]);
        Y=Y+ipermute(Y,[1,jj,ll,kk])+ipermute(Y,[1,jj,kk,ll]);
    end
    Y=reshape(Y,[nd,nz^3]);
end
end

function [Tz,others,eigval,retcode,options]=solve_first_order(structural_matrices,others,siz,pos,options,k_future)
dv=structural_matrices.dv;
Q=structural_matrices.transition_matrices.Q;

[dbf_plus,ds_0,dp_0,db_0,df_0,dpb_minus]=utils.solve.pull_first_order_partitions(dv,pos.v);

[Tz_pb,eigval,retcode,options]=dsge_solver_first_order_autoregress_h(dbf_plus,ds_0,dp_0,db_0,df_0,dpb_minus,Q,siz,pos,options);

Tz=cell(1,siz.h);
if ~retcode
    % (non-)certainty equivalence
    %----------------------------
    dt_t=zeros(siz.nd,siz.h);
    A0sig=zeros(siz.nd,siz.nT,siz.h);
    A0=zeros(siz.nd,siz.nT,siz.h);
    Tz_sig=zeros(siz.nT,siz.h);
    df_Lf_Tzp_Lp=zeros(siz.nd,siz.nT);
    db_Lb_Tzb_Lb=zeros(siz.nd,siz.nT);
    UUi=zeros(siz.nd,siz.nT,siz.h);
    others.dbf_plus=cell(siz.h);
    Tz_e_rt=cell(1,siz.h);
    for rt=1:siz.h
        de_0_rt=0;
        for rplus=1:siz.h
            ds_0=dv{rt,rplus}(:,pos.v.s_0);
            dp_0=dv{rt,rplus}(:,pos.v.p_0);
            db_0=dv{rt,rplus}(:,pos.v.b_0);
            df_0=dv{rt,rplus}(:,pos.v.f_0);
            A0_0_1=[ds_0,dp_0,db_0,df_0];
            A0(:,:,rt)=A0(:,:,rt)+A0_0_1;
            % provision for non-certainty equivalence
            %----------------------------------------
            dtheta_plus=dv{rt,rplus}(:,pos.v.theta_plus);
            df_plus=dv{rt,rplus}(:,pos.v.f_plus);
            db_plus=dv{rt,rplus}(:,pos.v.b_plus);
            df_Lf_Tzp_Lp(:,pos.t.p)=df_plus*Tz_pb(pos.t.f,1:siz.np,rplus);% place in the p position
            db_Lb_Tzb_Lb(:,pos.t.b)=db_plus*Tz_pb(pos.t.b,siz.np+(1:siz.nb),rplus);% place in the b position
            A0sig(:,:,rt) = A0sig(:,:,rt) + A0_0_1 + df_Lf_Tzp_Lp + db_Lb_Tzb_Lb;
            dt_t(:,rt)=dt_t(:,rt)+dtheta_plus*others.theta_hat{rt,rplus};
            
            % provision for shock impacts
            %----------------------------
            de_0=dv{rt,rplus}(:,pos.v.e_0);
            others.dbf_plus{rt,rplus}=dv{rt,rplus}(:,pos.v.bf_plus);
            UUi(:,:,rt)=UUi(:,:,rt)+A0_0_1;
            UUi(:,pos.t.pb,rt)=UUi(:,pos.t.pb,rt)+others.dbf_plus{rt,rplus}*Tz_pb(pos.t.bf,:,rplus);
            de_0_rt=de_0_rt+de_0;
        end
        % shock impacts (current)
        %------------------------
        UUi(:,:,rt)=UUi(:,:,rt)\eye(siz.nT);
        Tz_e_rt{rt}=-UUi(:,:,rt)*de_0_rt;
        Tz{rt}=[Tz_pb(:,:,rt),Tz_sig(:,rt),Tz_e_rt{rt}];
    end
    
    % shock impacts (future)
    %-----------------------
    for ik=1:k_future
        Tz_e_r0=Tz_e_rt;
        for rt=1:siz.h
            Sdbf_plus_rt_Tz_e0=0;
            for rplus=1:siz.h
                Sdbf_plus_rt_Tz_e0=Sdbf_plus_rt_Tz_e0+others.dbf_plus{rt,rplus}*Tz_e_r0{rplus}(pos.t.bf,:);
            end
            Tz_e_rt{rt}=-UUi(:,:,rt)*Sdbf_plus_rt_Tz_e0;
            Tz{rt}=[Tz{rt},Tz_e_rt{rt}];
        end
    end
    
    % now solve sum(A+*Tz_sig(+)+A0_sig*Tz_sig+dt_t=0
    %-------------------------------------------------
    % first we augment dt_t with the user_resid so that the first-order
    % approximation recoups the zero-th order approximation.
    dt_t=dt_t+structural_matrices.user_resids;
    Tz_sig=solve_perturbation_impact(Tz_sig,A0sig,others.dbf_plus,dt_t);
    if any(Tz_sig(:))
        for rt=1:siz.h
            Tz{rt}(:,siz.np+siz.nb+1)=Tz_sig(:,rt);
        end
    end
    others.Ui=UUi;
end

    function Tz_sig=solve_perturbation_impact(Tz_sig,A0sig,dbf_plus,dt_t)
        if any(dt_t(:))% then automatically h>1
            % use a qr decomposition to solve a small system. Given the structure
            % of the system, it is enough to precondition it.
            %-----------------------------------------------
            for r0=1:siz.h
                A0_sig_i=A0sig(:,:,r0)\eye(siz.nd);
                for r1=1:siz.h
                    dbf_plus{r0,r1}=A0_sig_i*dbf_plus{r0,r1};
                end
                dt_t(:,r0)=A0_sig_i*dt_t(:,r0);
            end
            % now we solve the system sum(A+*Tz_sig(+)+Tz_sig+dt_t=0 first for
            % variables p,b,f and then for variables s
            clear A0sig
            
            % solve the small system without static variables
            %------------------------------------------------
            % the direct solution implemented below is not efficient in very large
            % systems...
            npbf=siz.np+siz.nb+siz.nf;
            A=zeros(npbf*siz.h);
            for r0=1:siz.h
                row_=(r0-1)*npbf+1:r0*npbf;
                for r1=1:siz.h
                    col_=(r1-1)*npbf+1:r1*npbf;
                    A(row_,col_)=[zeros(npbf,siz.np),dbf_plus{r0,r1}(siz.ns+1:end,:)];
                    if r0==r1
                        A(row_,col_)=A(row_,col_)+eye(npbf);
                    end
                end
            end
            Tz_sig_PBF=uminus(dt_t(siz.ns+1:end,:));
            Tz_sig_PBF=A\Tz_sig_PBF(:);
            Tz_sig(siz.ns+1:end,:)=reshape(Tz_sig_PBF,npbf,siz.h);
            
            % solve the static variables
            %---------------------------
            for r0=1:siz.h
                Tz_sig(1:siz.ns,r0)=dt_t(1:siz.ns,r0);
                for r1=1:siz.h
                    Tz_sig(1:siz.ns,r0)=Tz_sig(1:siz.ns,r0)+...
                        dbf_plus{r0,r1}(1:siz.ns,:)*Tz_sig(pos.t.bf,r1);
                end
                Tz_sig(1:siz.ns,r0)=uminus(Tz_sig(1:siz.ns,r0));
            end
        end
    end
end

function [Tz_pb,eigval,retcode]=dsge_solver_first_order_autoregress_1(...
    dbf_plus,ds_0,dp_0,db_0,df_0,dpb_minus,siz,options)
nregs=numel(dbf_plus);
tolerance=1e-9;
all_same=true;
if nregs>1
    ncols_bf=siz.nb+siz.nf;
    ncols_pb=siz.np+siz.nb;
    if ncols_bf
        lead_=@(x)vec(dbf_plus{x});
        lead_1=lead_(1);
    end
    curr_=@(x)vec([ds_0{x},dp_0{x},db_0{x},df_0{x}]);
    curr_1=curr_(1);
    if ncols_pb
        lag_=@(x)vec(dpb_minus{x});
        lag_1=lag_(1);
    end
    for ireg=2:nregs
        if ncols_bf
            all_same=all_same && max(abs(lead_1-lead_(ireg)))<tolerance;
        end
        all_same=all_same && max(abs(curr_1-curr_(ireg)))<tolerance;
        if ncols_pb
            all_same=all_same && max(abs(lag_1-lag_(ireg)))<tolerance;
        end
        if ~all_same
            break
        end
    end
end

if ~all_same % a diagonal transition matrix with entirely different regimes
    % solve one at a time
    retcode=0;
    nvars=siz.ns+siz.np+siz.nb+siz.nf;
    Tz_pb=zeros(nvars,siz.np+siz.nb,nregs);
    eigval=cell(1,nregs);
    for ireg=1:nregs
        if ~retcode
            [Tsol,eigval{ireg},retcode]=dsge_solver_first_order_autoregress_1(...
                dbf_plus(ireg),ds_0(ireg),dp_0(ireg),db_0(ireg),df_0(ireg),...
                dpb_minus(ireg),siz,options);
            if ~retcode
                Tz_pb(:,:,ireg)=Tsol;
            end
        end
    end
    if ~retcode
        eigval=cell2mat(eigval);
    end
    return
end

rise_qz_criterium=sqrt(eps);
switch lower(options.solver)
    case {'rise_1'}
        [Tz_pb,eigval,retcode]=rise_solve_constant();
    case {'klein'}
        [Tz_pb,eigval,retcode]=klein_solve();
    case {'aim'}
        [Tz_pb,eigval,retcode]=aim_solve();
    case {'sims'}
        [Tz_pb,eigval,retcode]=sims_solve();
    otherwise
        error(['unknown solver ',parser.any2str(options.solver)])
end

if nregs>1
    Tz_pb=Tz_pb(:,:,ones(1,nregs));
end

    function [Tz_pb,eigval,retcode]=aim_solve(varargin)
        error('the aim solver is not yet implemented')
    end

    function [Tz_pb,eigval,retcode]=sims_solve(varargin)
        error('the sims solver is not yet implemented')
    end

    function [TT,SS,Z,eigval,retcode]=process_eigenvalues(TT,SS,Q,Z,npred)
        % Ordered inverse eigenvalues
        %----------------------------
        eigval = ordeig(TT,SS);
        stable = abs(eigval) >= 1 + rise_qz_criterium;
        nstable = sum(stable);
        unit = abs(abs(eigval)-1) < rise_qz_criterium;
        nunit = sum(unit);
        
        retcode=0;
        if nstable+nunit<npred
            retcode=22; % no solution
        elseif nstable+nunit>npred
            retcode=21; % multiple solutions
        else
            % Clusters of unit, stable, and unstable eigenvalues.
            clusters = zeros(size(eigval));
            
            % Unit roots first.
            %------------------
            clusters(unit) = 2;
            
            % Stable roots second.
            %---------------------
            clusters(stable) = 1;
            
            % Unstable roots last.
            %---------------------
            
            % Re-order by the clusters.
            %--------------------------
            [TT,SS,~,Z] = ordqz(TT,SS,Q,Z,clusters);
        end
        % Undo the eigval inversion.
        %---------------------------
        infeigval = eigval == 0;
        eigval(~infeigval) = 1./eigval(~infeigval);
        eigval(infeigval) = Inf;
    end

    function [Tz_pb,eigval,retcode]=klein_solve()
        % put system in the form a*x(t+1)=b*x(t) where x=[x0,xf];
        %--------------------------------------------------------
        nbf=siz.nb+siz.nf;
        bf_loc=siz.ns+siz.np+(1:siz.nb+siz.nf);
        pb_loc=siz.ns+(1:siz.np+siz.nb);
        B0=[ds_0{1,1},dp_0{1,1},db_0{1,1},df_0{1,1}];
        Bminus=[zeros(siz.nd,siz.ns),dpb_minus{1,1},zeros(siz.nd,siz.nf)];
        
        a=[B0,dbf_plus{1,1}
            zeros(nbf,siz.nd+nbf)];
        a(siz.nd+1:end,bf_loc)=eye(nbf);
        
        b=[Bminus,zeros(siz.nd,nbf)
            zeros(nbf,siz.nd),-eye(nbf)];
        b=-b;
        
        [Tz_pb,eigval,retcode] = solab(siz.nd);
        if ~retcode
            Tz_pb=Tz_pb(:,pb_loc);
        end
        
        function [sol,eigval,retcode] = solab(npred)
            % npred = number of stable guys
            [TT,SS,Q,Z] = qz(full(a),full(b));      % upper triangular factorization of the matrix pencil b-za
            
            % process eigenvalues
            %---------------------
            [TT,SS,Z,eigval,retcode]=process_eigenvalues(TT,SS,Q,Z,npred);
            
            sol=[];
            if ~retcode
                z11 = Z(1:npred,1:npred);
                
                z11i = z11\eye(npred);
                s11 = TT(1:npred,1:npred);
                t11 = SS(1:npred,1:npred);
                
                dyn = s11\t11;
                sol = real(z11*dyn*z11i);
                % z21 = Z(npred+1:end,1:npred);
                % f = real(z21*z11i); % already included in the lower part of p
            end
        end
    end

    function [Tzp,eigval,retcode]=rise_solve_constant()
        % state variables (lags): pred,both
        %----------------------------------
        Apb_minus=[dpb_minus{1,1}
            sparse(siz.nb,siz.nb+siz.np)]; % auxiliary equations for 'both' variables
        
        % forward-looking variables (leads): static,both,forward
        %-------------------------------------------------------
        Asbf_plus=[sparse(siz.nd,siz.ns),dbf_plus{1,1}
            sparse(siz.nb,siz.nb+siz.nf+siz.ns)]; % auxiliary equations for 'both' variables
        
        % forward-looking variables (current): static,both,forward
        %---------------------------------------------------------
        Asbf_0=[ds_0{1,1},sparse(siz.nd,siz.nb),df_0{1,1}
            sparse(siz.nb,siz.ns),speye(siz.nb),sparse(siz.nb,siz.nf)]; % auxiliary equations for 'both' variables
        
        % state variables (current): pred,both
        %-------------------------------------
        Apb_0=[dp_0{1,1},db_0{1,1}
            sparse(siz.nb,siz.np),-speye(siz.nb)]; % auxiliary equations for 'both' variables
        [Tzp,eigval,retcode]=rise_solve_1(Asbf_plus,Apb_0,Asbf_0,Apb_minus);
        if ~retcode
            % Re-order [s,b,f,p,b] as [s,p,b,f].we simply can ignore the last b
            %------------------------------------------------------------------
            static_=1:siz.ns;
            pred_=siz.ns+siz.nb+siz.nf+(1:siz.np);
            both_=siz.ns+(1:siz.nb);
            frwrd_=siz.ns+siz.nb+(1:siz.nf);
            order_var=[static_,pred_,both_,frwrd_];
            Tzp=Tzp(order_var,:);
        end
        function [Tzp,eigval,retcode]=rise_solve_1(Afrwrd_plus,Apred_0,Afrwrd_0,Apred_minus)
            A=[Apred_0,Afrwrd_plus]; % pred,frwrd
            B=-[Apred_minus,Afrwrd_0]; % pred,frwrd
            npred=size(Apred_0,2);
            nfrwrd=size(Afrwrd_0,2);
            % real schur decomposition
            %-------------------------
            [TT,SS,Q,Z] = qz(full(A),full(B),'real');
            % so we have Q*A*Z = TT, Q*B*Z = SS.
            
            % process eigenvalues
            %---------------------
            [TT,SS,Z,eigval,retcode]=process_eigenvalues(TT,SS,Q,Z,npred);
            
            Tzp=[];
            if ~retcode
                % define
                %-------
                W=Z.';
                % partition matrices
                %-------------------
                pred=1:npred;
                frwrd=npred+(1:nfrwrd);
                W11=W(pred,pred);
                W12=W(pred,frwrd);
                W21=W(frwrd,pred);
                W22=W(frwrd,frwrd);
                
                S11=SS(pred,pred); % S12=SS(pred,frwrd); % S21=SS(frwrd,pred); % S22=SS(frwrd,frwrd);
                
                T11=TT(pred,pred); % T12=TT(pred,frwrd); % T21=TT(frwrd,pred); % T22=TT(frwrd,frwrd);
                
                % form solution: forward-looking variables
                %-----------------------------------------
                Fzp=-W22\W21;
                tmp=W11+W12*Fzp;
                
                % form solution: predetermined variables
                %---------------------------------------
                Pzp=(T11*tmp)\S11*tmp;
                
                % final solution matrix: Forward-looking+predetermined
                %-----------------------------------------------------
                Tzp=[Fzp;Pzp];
                % check again
                if any(isnan(Tzp(:)))
                    retcode=22;
                end
            end
        end
    end
end

function [Tz_pb,eigval,retcode,options]=...
    dsge_solver_first_order_autoregress_h(dbf_plus,ds_0,dp_0,db_0,df_0,...
    dpb_minus,Q,siz,pos,options)

% options
%--------
bf_cols_adjusted=pos.t.bf;
pb_cols_adjusted=pos.t.pb;

nd_adjusted=siz.nd;
siz_adj=siz; % adjusted sizes
accelerate=options.solve_accelerate && siz.ns;
% aggregate A0 and A_
%--------------------
d0=cell(1,siz.h);
for r0=1:siz.h
    for r1=2:siz.h
        ds_0{r0,1}=ds_0{r0,1}+ds_0{r0,r1};
        dp_0{r0,1}=dp_0{r0,1}+dp_0{r0,r1};
        db_0{r0,1}=db_0{r0,1}+db_0{r0,r1};
        df_0{r0,1}=df_0{r0,1}+df_0{r0,r1};
        dpb_minus{r0,1}=dpb_minus{r0,1}+dpb_minus{r0,r1};
    end
    d0{r0}=[ds_0{r0,1},dp_0{r0,1},db_0{r0,1},df_0{r0,1}];
    % eliminate static variables for speed
    %-------------------------------------
    if accelerate
        if r0==1
            Abar_minus_s=cell(1,siz.h);
            R_s_s=cell(1,siz.h);
            R_s_ns=cell(1,siz.h);
            Abar_plus_s=cell(siz.h);
            nd_adjusted=nd_adjusted-siz.ns;
            bf_cols_adjusted=bf_cols_adjusted-siz.ns;
            pb_cols_adjusted=pb_cols_adjusted-siz.ns;
            siz_adj.ns=0;
            siz_adj.nd=siz_adj.nd-siz.ns;
            siz_adj.nT=siz_adj.nT-siz.ns;
        end
        [Q0,d0{r0}]=qr(d0{r0});
        dpb_minus{r0,1}=Q0'*dpb_minus{r0,1};
        Abar_minus_s{r0}=dpb_minus{r0,1}(1:siz.ns,:);
        dpb_minus{r0,1}=dpb_minus{r0,1}(siz.ns+1:end,:);
        for r1=1:siz.h
            dbf_plus{r0,r1}=Q0'*dbf_plus{r0,r1};
            Abar_plus_s{r0,r1}=dbf_plus{r0,r1}(1:siz.ns,:);
            dbf_plus{r0,r1}=dbf_plus{r0,r1}(siz.ns+1:end,:);
        end
        R_s_s{r0}=d0{r0}(1:siz.ns,1:siz.ns);
        R_s_ns{r0}=d0{r0}(1:siz.ns,siz.ns+1:end);
        d0{r0}=d0{r0}(siz.ns+1:end,siz.ns+1:end);
        ds_0{r0,1}=d0{r0}(:,1:siz_adj.ns);
        dp_0{r0,1}=d0{r0}(:,siz_adj.ns+(1:siz_adj.np));
        db_0{r0,1}=d0{r0}(:,siz_adj.ns+siz_adj.np+(1:siz_adj.nb));
        df_0{r0,1}=d0{r0}(:,siz_adj.ns+siz_adj.np+siz_adj.nb+(1:siz_adj.nf));
    end
end
ds_0=ds_0(:,1)';
dp_0=dp_0(:,1)';
db_0=db_0(:,1)';
df_0=df_0(:,1)';

if isempty(options.solver)
    is_evs=is_eigenvalue_solver();
    if is_evs
        options.solver='rise_1';
    else
        options.solver='mfi';
    end
end

T0=dsge_tools.utils.msre_initial_guess(d0,dpb_minus,dbf_plus,...
    options.solve_initialization);

model_class=isempty(pb_cols_adjusted)+2*isempty(bf_cols_adjusted);

retcode=0;
eigval=[];
switch model_class
    case 1 % forward-looking models
        Tz_pb=0*T0;
    case 2 % backward-looking models
        Tz_pb=dsge_tools.utils.msre_initial_guess(d0,dpb_minus,dbf_plus,'backward');
    case 3 % static models
        Tz_pb=0*T0;
    case 0 % hybrid models
        kron_method=strncmpi(options.solver,'mnk',3);
        
        if strcmpi(options.solver,'mfi')
            iterate_func=@(x)msre_solvers.functional_iteration_h(x,dbf_plus,d0,...
                dpb_minus,bf_cols_adjusted,pb_cols_adjusted);
        elseif any(strcmpi(options.solver,{'mnk','mn'}))
            iterate_func=@(x)msre_solvers.newton_iteration_h(x,dbf_plus,d0,dpb_minus,...
                bf_cols_adjusted,pb_cols_adjusted,kron_method,options);
        elseif any(strcmpi(options.solver,{'mfi_full','mnk_full','mn_full'}))
            [Gplus01,A0,Aminus,T0]=full_state_matrices(siz_adj,dbf_plus,d0,dpb_minus,T0);
            if strcmpi(options.solver,'mfi_full')
                iterate_func=@(x)msre_solvers.functional_iteration_h_full(x,Gplus01,A0,Aminus);
            else
                iterate_func=@(x)msre_solvers.newton_iteration_h_full(x,Gplus01,A0,Aminus,...
                    kron_method,options);
            end
        elseif strcmpi(options.solver,'fwz')
            [Gplus01,A0,Aminus,T0]=full_state_matrices(siz_adj,dbf_plus,d0,dpb_minus,T0);
            [iterate_func,solution_func,inverse_solution_func]= ...,sampling_func
                msre_solvers.fwz_newton_system(Gplus01,A0,Aminus,Q);
            T0=inverse_solution_func(T0);
        end
        
        switch lower(options.solver)
            case {'rise_1','klein','aim','sims'}
                dbf_plus_row=reconfigure_aplus();
                [Tz_pb,eigval,retcode]=dsge_solver_first_order_autoregress_1(...
                    dbf_plus_row,ds_0,dp_0,db_0,df_0,dpb_minus,siz_adj,options);
            case {'mfi','mfi_full','mnk','mnk_full','mn','mn_full','fwz'}
                [Tz_pb,~,retcode]=fix_point_iterator(iterate_func,T0,options);
                if  ~retcode && strcmpi(options.solver,'fwz')
                    Tz_pb=solution_func(Tz_pb);
                end
                if any(strcmpi(options.solver,{'mfi_full','mnk_full','mn_full'}))
                    Tz_pb=reshape(Tz_pb,[size(Tz_pb,1),size(Tz_pb,1),siz_adj.h]);
                end
                if any(strcmpi(options.solver,{'fwz','mfi_full','mnk_full','mn_full'}))
                    Tz_pb=Tz_pb(:,siz_adj.ns+(1:siz_adj.np+siz_adj.nb),:);
                end
            otherwise
                % user-defined solver
                %--------------------
                [Gplus01,A0,Aminus,T0]=full_state_matrices(siz_adj,dbf_plus,d0,dpb_minus,T0);
                
                [Tz_pb,~,retcode]=options.solver(Gplus01,A0,Aminus,Q,T0);
                
                % collect the relevant part
                %--------------------------
                if ~retcode
                    Tz_pb=Tz_pb(:,siz_adj.ns+(1:siz_adj.np+siz_adj.nb),:);
                end
        end
end


if ~retcode
    npb=siz.np+siz.nb;
    Tz_pb=reshape(Tz_pb,[nd_adjusted,npb,siz.h]);
    if accelerate
        % solve for the static variables
        %-------------------------------
        Sz_pb=zeros(siz.ns,npb,siz.h);
        for r0=1:siz.h
            ATT=0;
            for r1=1:siz.h
                ATT=ATT+Abar_plus_s{r0,r1}*Tz_pb(siz.np+1:end,:,r1); % <-- Tz_pb(npb+1:end,:,r1); we need also the both variables
            end
            ATT=ATT*Tz_pb(1:npb,:,r0);
            Sz_pb(:,:,r0)=-R_s_s{r0}\(Abar_minus_s{r0}+R_s_ns{r0}*Tz_pb(:,:,r0)+ATT);
        end
        Tz_pb=cat(1,Sz_pb,Tz_pb);
    end
end
    function Apl=reconfigure_aplus()
        Apl=cell(siz.h,1);
        for ii=1:siz.h
            if Q(ii,ii)
                Apl{ii}=dbf_plus{ii,ii}/Q(ii,ii);
            else
                error('knife-edge probability matrix prevents from recovering a matrix')
            end
        end
    end

    function flag=is_eigenvalue_solver()
        flag=true;
        if ~(siz.h==1||all(diag(Q)==1))
            d0_test=d0{1};
            dpb_minus_test=dpb_minus{1};
            dbf_plus0=reconfigure_aplus();
            dbf_plus_test=dbf_plus0{1};
            % check whether it is shocks only
            for st=2:siz.h
                t0=get_max(d0_test-d0{st});
                tplus=get_max(dbf_plus_test-dbf_plus0{st});
                tminus=get_max(dpb_minus_test-dpb_minus{st});
                tmax=max([t0,tplus,tminus]);
                if tmax>1e-9
                    flag=false;
                    break
                end
            end
        end
        function m=get_max(x)
            m=max(abs(x(:)));
        end
    end
end
