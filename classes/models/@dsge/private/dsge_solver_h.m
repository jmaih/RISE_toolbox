function [T,eigval,retcode,obj]=dsge_solver_h(obj,structural_matrices)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

% the obj going out probably contains the changed options

if isempty(obj)
    
    if nargout>1
        
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    
    end
    
    T=struct();
    
    return
    
end

%% begin
T=struct();

% initial condition for linear solvers
%-------------------------------------
x0_equal_D=false;

% pre-computation of the Czz, Czzz, Czzzz, Czzzzz terms
%------------------------------------------------------
precompute_czzzzz=false;

% options.solve_order 1
%--------
if obj.options.solve_order>=1
    
    debug=obj.options.debug;
    
    [pos,siz,shock_horizon]=dsge_tools.rehash_topology(obj,structural_matrices);
    
    % Structure of elements that will move across different orders
    %-------------------------------------------------------------
    others=struct();
    
    [T.Tz,others,eigval,retcode,obj.options]=solve_first_order(structural_matrices,...
        others,siz,pos,obj.options,shock_horizon);
    
    % higher orders
    %--------------
    if obj.options.solve_order>=2 && ~retcode
        
        if obj.options.solve_accelerate||debug
            
            [shrink,expand]=utils.kronecker.shrink_expand(siz.nz,oo);
            
        end
        
        % shortcuts to functions
        %-----------------------
        is_computable=@utils.cr.is_computable;
        A_times_kron_I_B=@utils.kronecker.A_times_kron_I_B;
        A_times_sum_perms=@utils.kronecker.A_times_sum_perms;
        Pfunc=@(ABCD,matsizes,varargin)...
            utils.kronecker.sum_permutations(ABCD,matsizes,...
            struct('use_old_algo',true),...
            varargin{:});
        kron_Q1_Qk_times_A=@utils.kronecker.kron_Q1_Qk_times_A;
        A_times_kron_Q1_Qk=@utils.kronecker.A_times_kron_Q1_Qk;
        A_times_k_kron_B=@utils.kronecker.A_times_k_kron_B;
        dv_vz_omega=@utils.cr.dv_vz_omega;
        kronall=@utils.kronecker.kronall;
        [T,retcode]=solve_higher_orders(T,others,obj.options.solve_accelerate);
    end
    
    % solve for growth constant
    %--------------------------
    if ~retcode
        T=growth_component_solver(obj,pos,T);
    end
end

    function [T,retcode]=solve_higher_orders(T,others,accelerate)
        
        % higher-order moments
        %----------------------
        [Eu{1:obj.options.solve_order}]=dsge_tools.u_higher_order_moments(siz);
        
        a0_z=sparse(siz.nv,siz.nz);
        a0_z(pos.v.b_minus,pos.z.b)=eye(siz.nb);
        a0_z(pos.v.p_minus,pos.z.p)=eye(siz.np);
        a0_z(pos.v.e_0,pos.z.e_0)=eye(siz.ne);
        a0_z=repmat({a0_z},siz.h,siz.h);
        a1_z=repmat({sparse(siz.nv,siz.nz)},1,siz.h);
        
        hz=sparse(siz.nz,siz.nz);
        hz(pos.z.sig,pos.z.sig)=1;
        hz(siz.np+siz.nb+1+(1:siz.ne*shock_horizon),pos.z.e_plus)=eye(shock_horizon*siz.ne);
        hz=repmat({hz},1,siz.h);
        
        dbf_plus=others.dbf_plus;
        for rt=1:siz.h
            for rplus=1:siz.h
                % this preconditioning could also be done once and for all and
                % should be included in Aplus
                dbf_plus{rt,rplus}=others.Ui(:,:,rt)*dbf_plus{rt,rplus};
            end
        end
        
%         tic
        Dzz=second_order_rhs();
        Czz=cell(1,siz.h);
        [T.Tzz,retcode]=solve_generalized_sylvester(Dzz,2);
        clear Dzz Czz
%         fprintf(1,'Order 2 done in %0.4f seconds \n\n',toc);
        
        if obj.options.solve_order>2 && ~retcode
            a0_zz=repmat({sparse(siz.nv,siz.nz^2)},siz.h,siz.h);
            a1_zz=repmat({sparse(siz.nv,siz.nz^2)},1,siz.h);
            hzz=repmat({sparse(siz.nz,siz.nz^2)},1,siz.h);
            
%             tic
            Dzzz=third_order_rhs();
            Czzz=cell(1,siz.h);
            [T.Tzzz,retcode]=solve_generalized_sylvester(Dzzz,3);
            clear Dzzz Czzz
%             fprintf(1,'Order 3 done in %0.4f seconds \n\n',toc);
            
            if obj.options.solve_order>3 && ~retcode
                a0_zzz=repmat({sparse(siz.nv,siz.nz^3)},siz.h,siz.h);
                a1_zzz=repmat({sparse(siz.nv,siz.nz^3)},1,siz.h);
                hzzz=repmat({sparse(siz.nz,siz.nz^3)},1,siz.h);
                Dzzzz=fourth_order_rhs();
                
%                 tic
                Czzzz=cell(1,siz.h);
                [T.Tzzzz,retcode]=solve_generalized_sylvester(Dzzzz,4);
                clear Dzzzz Czzzz
%                 fprintf(1,'Order 4 done in %0.4f seconds \n\n',toc);
                
                if obj.options.solve_order>4 && ~retcode
                    hzzzz=repmat({sparse(siz.nz,siz.nz^4)},1,siz.h);
                    a0_zzzz=sparse(siz.nv,siz.nz^4);
                    a1_zzzz=sparse(siz.nv,siz.nz^4);
                    Dzzzzz=fifth_order_rhs();
                    
%                     tic
                    Czzzzz=cell(1,siz.h);
                    [T.Tzzzzz,retcode]=solve_generalized_sylvester(Dzzzzz,5);
                    clear Dzzzzz
%                     fprintf(1,'Order 5 done in %0.4f seconds \n\n',toc);
                    
                    if obj.options.solve_order>5
                        error('perturbations of order greater than 5 not implemented');
                    end
                end
            end
        end
        
        function Dzz=second_order_rhs()
            Dzz=preallocate_rhs(2);
            for r0=1:siz.h
                hz{r0}(pos.z.pb,:)=T.Tz{r0}(pos.t.pb,:);
                for r1=1:siz.h
                    if r0==1
                        a1_z{r1}(pos.v.bf_plus,:)=T.Tz{r1}(pos.t.bf,:);
                    end
                    a0_z{r0,r1}(pos.v.t_0,:)=T.Tz{r0};
                    a0_z{r0,r1}(pos.v.bf_plus,:)=T.Tz{r1}(pos.t.bf,:)*hz{r0};
                    Dzz(:,:,r0)=Dzz(:,:,r0)+dvv_Evz_vz();
                end
                % precondition
                Dzz(:,:,r0)=-others.Ui(:,:,r0)*Dzz(:,:,r0);
            end
            
            function res=dvv_Evz_vz()
                res=sparse(siz.nd,siz.nz^2);
                if is_computable(structural_matrices.dvv{r0,r1})
                    res=res+A_times_k_kron_B(structural_matrices.dvv{r0,r1},a0_z{r0,r1},2);
                    if is_computable(a1_z{r1})
                        res=res+A_times_k_kron_B(structural_matrices.dvv{r0,r1},a1_z{r1},2)*Eu{2};
                    end
                end
            end
        end
        
        function Dzzz=third_order_rhs()
            Dzzz=preallocate_rhs(3);
            for r0=1:siz.h
                hzz{r0}(pos.z.pb,:)=T.Tzz{r0}(pos.t.pb,:);
                for r1=1:siz.h
                    % second order
                    a0_zz{r0,r1}(pos.v.t_0,:)=T.Tzz{r0};
                    a0_zz{r0,r1}(pos.v.bf_plus,:)=T.Tz{r1}(pos.t.bf,:)*hzz{r0}+T.Tzz{r1}(pos.t.bf,:)*kron(hz{r0},hz{r0});
                    if r0==1
                        a1_zz{r1}(pos.v.bf_plus,:)=T.Tzz{r1}(pos.t.bf,:);
                    end
                    
                    Dzzz(:,:,r0)=Dzzz(:,:,r0)+dvvv_Evz_vz_vz();
                    
                    Dzzz(:,:,r0)=Dzzz(:,:,r0)+dvv_Evz_vzz();
                    
                    Dzzz(:,:,r0)=Dzzz(:,:,r0)+others.dbf_plus{rt,rplus}*Tzz_hz_hzz();
                end
                % precondition
                Dzzz(:,:,r0)=-others.Ui(:,:,r0)*Dzzz(:,:,r0);
            end
            
            function res=dvv_Evz_vzz()
                res=A_times_kron_Q1_Qk(...
                    structural_matrices.dvv{r0,r1},...
                    a0_z{r0,r1},a0_zz{r0,r1}+a1_zz{r1}*Eu{2});
                
                res=res+A_times_sum_perms(...
                    A_times_kron_Q1_Qk(structural_matrices.dvv{r0,r1},...
                    a1_z{r1},a1_zz{r1}),...
                    {Eu{2},hz{r0}},...
                    siz.nz*ones(3,2),...
                    true,[1,3,2]);
                
                res=dv_vz_omega(res,siz.nz,1);
            end
            
            function res=Tzz_hz_hzz()
                if is_computable(T.Tzz{r1}(pos.t.bf,:),hz{r0},hzz{r0})
                    res=A_times_kron_Q1_Qk(T.Tzz{r1}(pos.t.bf,:),hz{r0},hzz{r0});
                    res=dv_vz_omega(res,siz.nz,1);
                    % res=fvv_vx_vxx_omega_1(T.Tzz{r1}(pos.t.bf,:),hz{r0},hzz{r0});
                else
                    res=sparse(siz.nb+siz.nf,siz.nz^3);
                end
            end
            
            function res=dvvv_Evz_vz_vz()
                res=fvvv_vx_vx_vx(structural_matrices.dvvv{r0,r1},a0_z{r0,r1});
                
                res=res+A_times_sum_perms(...
                    structural_matrices.dvvv{r0,r1},...
                    {kron_Q1_Qk_times_A(Eu{2},a1_z{r1},a1_z{r1}),a0_z{r0,r1}},...
                    ones(3,1)*[siz.nv,siz.nz],...
                    true,[1,3,2],[3,1,2]);
            end
        end
        
        function Dzzzz=fourth_order_rhs()
            Dzzzz=preallocate_rhs(4);
            for r0=1:siz.h
                hzzz{r0}(pos.z.pb,:)=T.Tzzz{r0}(pos.t.pb,:);
                for r1=1:siz.h
                    % third order
                    if r0==1
                        a1_zzz{r1}(pos.v.bf_plus,:)=T.Tzzz{r1}(pos.t.bf,:);
                    end
                    a0_zzz{r0,r1}(pos.v.t_0,:)=T.Tzzz{r0};
                    a0_zzz{r0,r1}(pos.v.bf_plus,:)=...
                        fvvv_vx_vx_vx(T.Tzzz{r1}(pos.t.bf,:),hz{r0})+...
                        fvv_vx_vxx_omega_1(T.Tzz{r1}(pos.t.bf,:),hz{r0},hzz{r0})+...
                        T.Tz{r1}(pos.t.bf,:)*hzzz{r0};
                    
                    Dzzzz(:,:,r0)=Dzzzz(:,:,r0)+dvvvv_Evz_vz_vz_vz();
                    
                    Dzzzz(:,:,r0)=Dzzzz(:,:,r0)+dvvv_Evz_vz_vzz();
                    
                    Dzzzz(:,:,r0)=Dzzzz(:,:,r0)+dvv_Evz_vzzz();
                    
                    Dzzzz(:,:,r0)=Dzzzz(:,:,r0)+dvv_Evzz_vzz();
                    
                    Dzzzz(:,:,r0)=Dzzzz(:,:,r0)+others.dbf_plus{rt,rplus}*lambda_bf_XI01_4();
                end
                % precondition
                Dzzzz(:,:,r0)=-others.Ui(:,:,r0)*Dzzzz(:,:,r0);
            end
            
            function res=dvvvv_Evz_vz_vz_vz()
                res=A_times_k_kron_B(structural_matrices.dvvvv{r0,r1},a0_z{r0,r1},4);
                
                res=res+A_times_k_kron_B(structural_matrices.dvvvv{r0,r1},a1_z{r1},4)*Eu{4};
                
                tmp=kron_Q1_Qk_times_A(Eu{2},a1_z{r1},a1_z{r1});
                matsizes=ones(4,1)*[siz.nv,siz.nz];
                res=res+A_times_sum_perms(...
                    structural_matrices.dvvvv{r0,r1},...
                    {a0_z{r0,r1},a0_z{r0,r1},tmp},...
                    matsizes,true,...
                    [1,3,2,4],[1,3,4,2],[3,1,4,2],[3,4,1,2],[3,1,2,4]);
            end
            
            function res=dvvv_Evz_vz_vzz()
                res=A_times_kron_Q1_Qk(structural_matrices.dvvv{r0,r1},...
                    a0_z{r0,r1},a0_z{r0,r1},a0_zz{r0,r1});
                
                res=res+A_times_kron_Q1_Qk(structural_matrices.dvvv{r0,r1},...
                    a0_z{r0,r1},a0_z{r0,r1},a1_zz{r1}*Eu{2});
                
                res=res+A_times_kron_Q1_Qk(structural_matrices.dvvv{r0,r1},...
                    a1_z{r1},a1_z{r1},a1_zz{r1})*Eu{4};
                
                res=res+A_times_kron_Q1_Qk(structural_matrices.dvvv{r0,r1},...
                    kron_Q1_Qk_times_A(Eu{2},a1_z{r1},a1_z{r1}),a0_zz{r0,r1});
                
                matsizes=[ones(2,1)*[siz.nv,siz.nz];siz.nz(ones(2))];
                res0=A_times_sum_perms(...
                    structural_matrices.dvvv{r0,r1}*kron(speye(siz.nv^2),a1_zz{r1}),... A_times_kron_I_B(structural_matrices.dvvv{r0,r1},a1_zz{r1},siz.nv^2)
                    {a0_z{r0,r1},...
                    kron(a1_z{r1},speye(siz.nz))*Eu{2},hz{r0}},... kron_A_I_times_B(a1_z{r1},Eu{2},siz.nz)
                    matsizes,true,...
                    [1,2,4,3],[2,1,3,4],[2,1,4,3]);
                res=res+res0;
                res=dv_vz_omega(res,siz.nz,2);
                
                if debug
                    B=Pfunc(kron(Eu{2},hz{r0}),siz.nz(ones(3,2)),[1,3,2]);
                    A=kron_Q1_Qk_times_A(B,a1_z{r1},a1_zz{r1});
                    A=Pfunc(kron(a0_z{r0,r1},A),...
                        [ones(2,1)*[siz.nv,siz.nz];[siz.nv,siz.nz^2]],[2,1,3]);
                    res00=structural_matrices.dvvv{r0,r1}*A;
                    max(abs(res0(:)-res00(:)))
                    keyboard
                end
            end
            
            function res=dvv_Evz_vzzz()
                res=A_times_kron_Q1_Qk(structural_matrices.dvv{r0,r1},a0_z{r0,r1},a0_zzz{r0,r1});
                
                res=res+A_times_kron_Q1_Qk(structural_matrices.dvv{r0,r1},a1_z{r1},a1_zzz{r1})*Eu{4};
                
                res=res+A_times_sum_perms(...
                    A_times_kron_Q1_Qk(structural_matrices.dvv{r0,r1},a1_z{r1},a1_zzz{r1}),...
                    {Eu{2},hz{r0},hz{r0}},...
                    siz.nz(ones(4,2)),true,...
                    [1,3,2,4],[1,3,4,2]);
                
                res=res+A_times_kron_Q1_Qk(structural_matrices.dvv{r0,r1},...
                    a0_z{r0,r1},A_times_sum_perms(...
                    a1_zzz{r1},{hz{r0},Eu{2}},...
                    siz.nz(ones(3,2)),true,...
                    [2,1,3],[2,3,1]));
                
                omega1=utils.cr.omega(siz.nz,1);
                
                res0=A_times_kron_I_B(...
                    A_times_kron_Q1_Qk(...
                    A_times_kron_Q1_Qk(structural_matrices.dvv{r0,r1},...
                    a1_z{r1},a1_zz{r1}),...
                    Eu{2},hzz{r0}),...
                    omega1,siz.nz);
                res=res+res0;
                res=dv_vz_omega(res,siz.nz,3);
                
                if debug
                    res00=A_times_kron_Q1_Qk(...
                        A_times_kron_Q1_Qk(...
                        A_times_kron_Q1_Qk(structural_matrices.dvv{r0,r1},...
                        a1_z{r1},a1_zz{r1}),...
                        Eu{2},hzz{r0}),...
                        speye(siz.nz),...
                        omega1);
                    max(abs(res0(:)-res00(:)))
                    keyboard
                end
            end
            
            function res=dvv_Evzz_vzz()
                res=A_times_kron_Q1_Qk(structural_matrices.dvv{r0,r1},a0_zz{r0,r1},a0_zz{r0,r1});
                
                res=res+A_times_sum_perms(...
                    structural_matrices.dvv{r0,r1},...
                    {a0_zz{r0,r1},a1_zz{r1}*Eu{2}},...
                    ones(2,1)*[siz.nv,siz.nz^2],true,...
                    [2,1]);
                
                res=res+A_times_kron_Q1_Qk(...
                    structural_matrices.dvv{r0,r1},a1_zz{r1},a1_zz{r1})*Eu{4};
                
                res=res+A_times_sum_perms(...
                    A_times_kron_Q1_Qk(...
                    structural_matrices.dvv{r0,r1},a1_zz{r1},a1_zz{r1}),...
                    {hz{r0},Eu{2},hz{r0}},...
                    siz.nz(ones(4,2)),...
                    true,[1,2,4,3],...
                    [2,1,3,4],[2,1,4,3]);
                
                res=dv_vz_omega(res,siz.nz,4);
            end
            
            function res=lambda_bf_XI01_4()
                res0=A_times_kron_Q1_Qk(T.Tzzz{r1}(pos.t.bf,:),hz{r0},hz{r0},hzz{r0});
                
                res0=res0+A_times_kron_Q1_Qk(T.Tzzz{r1}(pos.t.bf,:),Eu{2},hzz{r0});
                
                res=dv_vz_omega(res0,siz.nz,2);
                if debug
                    % the multiplication seems to introduce some
                    % substantial but still negligible noise
                    res__=res0*utils.cr.omega(siz.nz,2);
                    disp(max(max(abs(res__-res)))) %e.g. 5.9605e-08
                    keyboard
                end
                
                res=res+dv_vz_omega(A_times_kron_Q1_Qk(T.Tzz{r1}(pos.t.bf,:),...
                    hz{r0},hzzz{r0}),siz.nz,3);
                
                res=res+dv_vz_omega(...
                    A_times_kron_Q1_Qk(T.Tzz{r1}(pos.t.bf,:),hzz{r0},hzz{r0}),...
                    siz.nz,4);
            end
        end
        
        function Dzzzzz=fifth_order_rhs()
            Dzzzzz=preallocate_rhs(5);
            for r0=1:siz.h
                hzzzz{r0}(pos.z.pb,:)=T.Tzzzz{r0}(pos.t.pb,:);
                for r1=1:siz.h
                    % fourth order
                    a0_zzzz(pos.v.t_0,:)=T.Tzzzz{r0};
                    a0_zzzz(pos.v.bf_plus,:)=redo_a0_zzzz();
                    a1_zzzz(pos.v.bf_plus,:)=T.Tzzzz{r1}(pos.t.bf,:);
                    
                    Dzzzzz(:,:,r0)=Dzzzzz(:,:,r0)+dvvvvv_Evz_vz_vz_vz_vz();
                    
                    Dzzzzz(:,:,r0)=Dzzzzz(:,:,r0)+dvvvv_Evz_vz_vz_vzz();
                    
                    Dzzzzz(:,:,r0)=Dzzzzz(:,:,r0)+dvvv_Evz_vz_vzzz();
                    
                    Dzzzzz(:,:,r0)=Dzzzzz(:,:,r0)+dvvv_Evz_vzz_vzz();
                    
                    Dzzzzz(:,:,r0)=Dzzzzz(:,:,r0)+dvv_Evz_vzzzz();
                    
                    Dzzzzz(:,:,r0)=Dzzzzz(:,:,r0)+dvv_Evzz_vzzz();
                    
                    Dzzzzz(:,:,r0)=Dzzzzz(:,:,r0)+others.dbf_plus{rt,rplus}*lambda_bf_XI01();
                end
                % precondition
                Dzzzzz(:,:,r0)=-others.Ui(:,:,r0)*Dzzzzz(:,:,r0);
            end
            
            function res=redo_a0_zzzz()
                % builds on lambda_bf_XI01_4 but it is different!!!!
                %---------------------------------------------------
                res=A_times_kron_Q1_Qk(T.Tzzz{r1}(pos.t.bf,:),hz{r0},hz{r0},hzz{r0});
                res=res+A_times_kron_Q1_Qk(T.Tzzz{r1}(pos.t.bf,:),hz{r0},hz{r0},hzz{r0});
                res=dv_vz_omega(res,siz.nz,2);
                res=res+dv_vz_omega(A_times_kron_Q1_Qk(T.Tzz{r1}(pos.t.bf,:),...
                    hz{r0},hzzz{r0}),siz.nz,3);
                res=res+dv_vz_omega(A_times_kron_Q1_Qk(T.Tzz{r1}(pos.t.bf,:),hzz{r0},hzz{r0}),...
                    siz.nz,4);
                % add the others
                %---------------
                res=res+T.Tz{r1}(pos.t.bf,:)*hzzzz{r0};
                res=res+A_times_kron_Q1_Qk(T.Tzzzz{r1}(pos.t.bf,:),hz{r0},hz{r0},hz{r0},hz{r0});
            end
            
            function res=dvvvvv_Evz_vz_vz_vz_vz()
                res=A_times_kron_Q1_Qk(...
                    structural_matrices.dvvvvv{r0,r1},a0_z{r0,r1},a0_z{r0,r1},a0_z{r0,r1},a0_z{r0,r1},a0_z{r0,r1});
                
                tmp=kron_Q1_Qk_times_A(Eu{2},a1_z{r1},a1_z{r1});
                
                matsizes=ones(5,1)*[siz.nv,siz.nz];
                res=res+A_times_sum_perms(structural_matrices.dvvvvv{r0,r1},...
                    {a0_z{r0,r1},a0_z{r0,r1},a0_z{r0,r1},tmp},...
                    matsizes,...
                    true,...
                    [1,2,4,3,5],[1,2,4,5,3],[1,4,2,5,3],...
                    [1,4,5,2,3],[4,1,5,2,3],[4,5,1,2,3],[4,1,2,5,3],...
                    [4,1,2,3,5],[1,4,2,3,5]);
                
                % exploit sparsity by using kron directly
                %-----------------------------------------
                tmp=kronall(a1_z{r1},a1_z{r1},a1_z{r1},a1_z{r1})*Eu{4};%<--tmp=kron_Q1_Qk_times_A(Eu{4},a1_z{r1},a1_z{r1},a1_z{r1},a1_z{r1});
                res=res+A_times_sum_perms(structural_matrices.dvvvvv{r0,r1},...
                    {a0_z{r0,r1},tmp},...
                    matsizes,...
                    true,[2,1,3,4,5],[2,3,1,4,5],[2,3,4,1,5],[2,3,4,5,1]);
            end
            
            function res=dvvvv_Evz_vz_vz_vzz()
                res=A_times_kron_Q1_Qk(...
                    structural_matrices.dvvvv{r0,r1},a0_z{r0,r1},a0_z{r0,r1},a0_z{r0,r1},a0_zz{r0,r1});
                
                res=res+A_times_kron_Q1_Qk(...
                    structural_matrices.dvvvv{r0,r1},a0_z{r0,r1},a0_z{r0,r1},a0_z{r0,r1},a1_zz{r1}*Eu{2});
                
                tmp=A_times_kron_Q1_Qk(...
                    structural_matrices.dvvvv{r0,r1},a1_z{r1},a1_z{r1},a1_z{r1},a1_zz{r1});
                res=res+A_times_sum_perms(tmp,{Eu{4},hz{r0}},...
                    siz.nz*ones(5,2),...
                    true,[1,2,3,5,4]);
                
                res=res+A_times_sum_perms(...
                    structural_matrices.dvvvv{r0,r1},...
                    {a0_z{r0,r1},kron_Q1_Qk_times_A(Eu{2},a1_z{r1},a1_z{r1}),a0_zz{r0,r1}},...
                    [ones(3,1)*[siz.nv,siz.nz];[siz.nv,siz.nz^2]],...
                    true,[2,1,3,4],[2,3,1,4]);
                
                tmp=A_times_sum_perms(...
                    kron(a1_z{r1},a1_zz{r1}),{Eu{2},hz{r0}},...
                    siz.nz*ones(3,2),...
                    true,...
                    [1,3,2]);
                res=res+A_times_sum_perms(...
                    structural_matrices.dvvvv{r0,r1},...
                    {a0_z{r0,r1},a0_z{r0,r1},tmp},...
                    [ones(3,1)*[siz.nv,siz.nz];[siz.nv,siz.nz^2]],...
                    true,[1,3,2,4],[3,1,2,4]...
                    );
                
                % invoke sparsity
                %------------------
                tmp=kronall(a1_z{r1},a1_z{r1},a1_zz{r1})*Eu{4};%<--tmp=kron_Q1_Qk_times_A(Eu{4},a1_z{r1},a1_z{r1},a1_zz{r1});
                res=res+A_times_sum_perms(...
                    structural_matrices.dvvvv{r0,r1},...
                    {a0_z{r0,r1},tmp},...
                    [ones(3,1)*[siz.nv,siz.nz];[siz.nv,siz.nz^2]],...
                    true,[2,1,3,4],[2,3,1,4]);
                
                res=dv_vz_omega(res,siz.nz,5);
            end
            
            function res=dvvv_Evz_vz_vzzz()
                res=A_times_kron_Q1_Qk(...
                    structural_matrices.dvvv{r0,r1},a0_z{r0,r1},a0_z{r0,r1},...
                    a0_zzz{r0,r1}+...
                    A_times_sum_perms(a1_zzz{r1},{hz{r0},Eu{2}},siz.nz*ones(3,2),...
                    true,[2,1,3],[2,3,1])...
                    );
                
                res=res+A_times_kron_Q1_Qk(...
                    structural_matrices.dvvv{r0,r1},...
                    kron_Q1_Qk_times_A(Eu{2},a1_z{r1},a1_z{r1}),a0_zzz{r0,r1});
                
                tmp=A_times_kron_Q1_Qk(...
                    structural_matrices.dvvv{r0,r1},...
                    a1_z{r1},a1_z{r1},a1_zzz{r1});
                res=res+A_times_sum_perms(...
                    tmp,{Eu{4},hz{r0}},...
                    [[siz.nz^2,siz.nz^2];siz.nz*ones(3,2)],true,...
                    [1,2,4,3],[1,4,2,3]);
                
                tmp=A_times_sum_perms(kron(a1_z{r1},a1_zzz{r1}),...
                    {Eu{2},hz{r0},hz{r0}},...
                    siz.nz*ones(4,2),...
                    true,[1,3,2,4],[1,3,4,2]);
                res=res+A_times_sum_perms(structural_matrices.dvvv{r0,r1},...
                    {a0_z{r0,r1},tmp},...
                    [ones(2,1)*[siz.nv,siz.nz];[siz.nv,siz.nz^3]],...
                    true,[2,1,3]);
                
                res=res+A_times_sum_perms(structural_matrices.dvvv{r0,r1},...
                    {a0_z{r0,r1},kron_Q1_Qk_times_A(Eu{4},a1_z{r1},a1_zzz{r1})},...
                    [ones(2,1)*[siz.nv,siz.nz];[siz.nv,siz.nz^3]],...
                    true,[2,1,3]);
                
                omega1=utils.cr.omega(siz.nz,1);
                tmp=kron_Q1_Qk_times_A(kron(Eu{2},hzz{r0}),a1_z{r1},a1_zz{r1});
                tmp=tmp*kron(speye(siz.nz),omega1);
                res=res+A_times_sum_perms(...
                    structural_matrices.dvvv{r0,r1},...
                    {a0_z{r0,r1},tmp},...
                    [ones(2,1)*[siz.nv,siz.nz];[siz.nv,siz.nz^3]],...
                    true,[2,1,3]);
                
                res=dv_vz_omega(res,siz.nz,6);
            end
            
            function res=dvvv_Evz_vzz_vzz()
                res=A_times_kron_Q1_Qk(...
                    structural_matrices.dvvv{r0,r1},a0_z{r0,r1},...
                    a0_zz{r0,r1},a0_zz{r0,r1}+a1_zz{r1}*Eu{2});
                
                tmp=A_times_sum_perms(kron(a1_zz{r1},a1_zz{r1}),...
                    {hz{r0},Eu{2},hz{r0}},...
                    siz.nz*ones(4,2),true,...
                    [1,2,4,3],[2,1,4,3],[2,1,3,4]);
                res=res+A_times_kron_Q1_Qk(...
                    structural_matrices.dvvv{r0,r1},a0_z{r0,r1},tmp); clear tmp
                
                res=res+A_times_kron_Q1_Qk(...
                    structural_matrices.dvvv{r0,r1},a0_z{r0,r1},...
                    kron(a1_zz{r1}*Eu{2},a0_zz{r0,r1})+...
                    kron(a1_zz{r1},a1_zz{r1})*Eu{4});
                
                res=res+A_times_sum_perms(structural_matrices.dvvv{r0,r1},...
                    {kron(a1_z{r1},a1_zz{r1})*kron(Eu{2},hz{r0}),a0_zz{r0,r1}},...
                    [[siz.nv,siz.nz];ones(2,1)*[siz.nv,siz.nz^2]],...
                    true,[1,3,2]);
                
                tmp=A_times_sum_perms(kronall(a1_z{r1},a1_zz{r1},a1_zz{r1}),...
                    {hz{r0},Eu{4}},siz.nz*ones(5,2),false,...
                    [2,3,1,4,5],[2,1,3,4,5]);
                res=res+A_times_sum_perms(structural_matrices.dvvv{r0,r1},...
                    {tmp},[[siz.nv,siz.nz];ones(2,1)*[siz.nv,siz.nz^2]],...
                    true,[1,3,2]);
                
                res=dv_vz_omega(res,siz.nz,7);
            end
            
            function res=dvv_Evz_vzzzz()
                res=A_times_kron_Q1_Qk(...
                    structural_matrices.dvv{r0,r1},a0_z{r0,r1},a0_zzzz+...
                    a1_zzzz*(Eu{4}+...
                    Pfunc(kronall(Eu{2},hz{r0},hz{r0}),siz.nz*ones(4,2),...
                    [1,3,2,4],[1,3,4,2],[3,1,4,2],[3,4,1,2],[3,1,2,4]))+...
                    dv_vz_omega(A_times_kron_Q1_Qk(a1_zzz{r1},Eu{2},hzz{r0}),siz.nz,2));
                
                res=res+A_times_kron_Q1_Qk(...
                    structural_matrices.dvv{r0,r1},a1_z{r1},a1_zzzz)*(...
                    Pfunc(kronall(Eu{2},hz{r0},hz{r0},hz{r0}),siz.nz*ones(5,2),...
                    [1,3,2,4,5],[1,3,4,2,5],[1,3,4,5,2])+...
                    Pfunc(kron(Eu{4},hz{r0}),siz.nz*ones(5,2),...
                    [1,2,3,5,4],[1,2,5,3,4],[1,5,2,3,4])...
                    );
                
                omega5=utils.cr.omega(siz.nz,5);
                tmp=A_times_kron_Q1_Qk(...
                    structural_matrices.dvv{r0,r1},a1_z{r1},a1_zzz{r1});
                res=res+A_times_sum_perms(tmp,{kronall(Eu{2},hz{r0},hzz{r0})*omega5},...
                    [siz.nz*ones(3,2);[siz.nz,siz.nz^2]],...
                    true,[1,3,2,4]); clear tmp
                
                omega6=utils.cr.omega(siz.nz,6);
                tmp=A_times_kron_Q1_Qk(structural_matrices.dvv{r0,r1},a1_z{r1},a1_zz{r1});
                res=res+A_times_kron_Q1_Qk(tmp,Eu{2},hzzz{r0})*omega6; clear tmp
                
                res=dv_vz_omega(res,siz.nz,8);
            end
            
            function res=dvv_Evzz_vzzz()
                res=A_times_kron_Q1_Qk(...
                    structural_matrices.dvv{r0,r1},a0_zz{r0,r1},...
                    a0_zzz{r0,r1}+...
                    a1_zzz{r1}*Pfunc(kron(hz{r0},Eu{2}),...
                    siz.nz*ones(3,2),[2,1,3],[2,3,1])...
                    );
                
                res=res+A_times_kron_Q1_Qk(...
                    structural_matrices.dvv{r0,r1},a1_zz{r1}*Eu{2},a0_zzz{r0,r1});
                
                res=res+A_times_kron_Q1_Qk(...
                    structural_matrices.dvv{r0,r1},a1_zz{r1},a1_zzz{r1})*(...
                    Pfunc(kronall(hz{r0},Eu{2},hz{r0},hz{r0}),siz.nz*ones(5,2),...
                    [1,2,4,3,5],[1,2,4,5,3],[2,1,3,4,5],[2,1,4,3,5],...
                    [2,1,4,5,3])+...
                    Pfunc(kron(hz{r0},Eu{4}),...
                    [siz.nz*ones(2,2);[siz.nz^3,siz.nz^3]],...
                    [2,1,3])+...
                    Pfunc(kron(Eu{4},hz{r0}),...
                    [[siz.nz^2,siz.nz^2];siz.nz*ones(3,2)],...
                    [1,2,4,3],[1,4,2,3]));
                
                omega5=utils.cr.omega(siz.nz,5);
                tmp=A_times_kron_Q1_Qk(...
                    structural_matrices.dvv{r0,r1},a1_zz{r1},a1_zz{r1});
                res=res+A_times_sum_perms(tmp,...
                    {kronall(hz{r0},Eu{2},hzz{r0})*omega5},...
                    [siz.nz*ones(3,2);[siz.nz,siz.nz^2]],...
                    true,[2,1,3,4]);
                
                res=dv_vz_omega(res,siz.nz,9);
            end
            
            function res=lambda_bf_XI01()
                res=dv_vz_omega(A_times_kron_Q1_Qk(...
                    T.Tzzzz{r1}(pos.t.bf,:),hz{r0},hz{r0},hz{r0},hzz{r0}),...
                    siz.nz,5);
                
                omega5=utils.cr.omega(siz.nz,5);
                res=res+A_times_sum_perms(...
                    T.Tzzzz{r1}(pos.t.bf,:),{Eu{2},hz{r0},hzz{r0}},...
                    [siz.nz*ones(3,2);[siz.nz,siz.nz^2]],...
                    true,[1,3,2,4],[3,1,2,4])*omega5;
                
                res=res+dv_vz_omega(A_times_kron_Q1_Qk(...
                    T.Tzzz{r1}(pos.t.bf,:),hz{r0},hz{r0},hzzz{r0}),siz.nz,6);
                
                res=res+dv_vz_omega(A_times_kron_Q1_Qk(...
                    T.Tzzz{r1}(pos.t.bf,:),Eu{2},hzzz{r0}),siz.nz,6);
                
                res=res+dv_vz_omega(A_times_kron_Q1_Qk(...
                    T.Tzzz{r1}(pos.t.bf,:),hz{r0},hzz{r0},hzz{r0}),siz.nz,7);
                
                res=res+dv_vz_omega(T.Tzz{r1}(pos.t.bf,:)*...
                    kron(hz{r0},hzzzz{r0}),siz.nz,8);
                
                res=res+dv_vz_omega(T.Tzz{r1}(pos.t.bf,:)*...
                    kron(hzz{r0},hzzz{r0}),siz.nz,9);
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

                nkept=sum(shrink{oo});
                if debug
                    Dtest=D;
                end
                D=D(:,shrink{oo},:);
                if debug
                    Dbig=D(:,expand{oo},:);
                    max(abs(Dtest(:)-Dbig(:)))
                end
            else
                if debug

                    D=D(:,shrink{oo},:);
                    D=D(:,expand{oo},:);
                end
                nkept=siz.nz^oo;
            end
            
            % expand final result: maybe, maybe not
            %--------------------------------------
            if x0_equal_D
                x0=D(:);
            else
                x0=[];
            end
            
            [X,retcode]=utils.optim.linear_systems_solver(@afun,D(:),x0,obj.options);
                
            if ~retcode
                if accelerate
                    X=reshape(X,[siz.nd,nkept,siz.h]);
                    X=X(:,expand{oo},:);
                else
                    X=reshape(X,[siz.nd,siz.nz^oo,siz.h]);
                end
                tmp=cell(1,siz.h);
                for r0=1:siz.h
                    tmp{r0}=sparse(X(:,:,r0));
                end
                X=tmp;
            end
            
            function ATC_plus_T=afun(tau)
                tau=reshape(tau,[siz.nd,nkept,siz.h]);
                ATC_plus_T=zeros(size(tau));
                for r00=1:siz.h
                    AT=0;
                    for r1=1:siz.h
                        AT=AT+dbf_plus{r00,r1}*tau(pos.t.bf,:,r1);
                    end
                    if accelerate
                        AT=AT(:,expand{oo});
                    end
                    ATCzzz=temporary_term();
                    ATC_plus_T(:,:,r00)=ATCzzz+tau(:,:,r00);
                end
                ATC_plus_T=ATC_plus_T(:);
                function ATC=temporary_term()
                    % core terms
                    %-----------
                    AT=sparse(AT);
                    if oo==2
                        if precompute_czzzzz
                            if isempty(Czz{r00})
                                Czz{r00}=kron(hz{r00},hz{r00})+Eu{oo};
                            end
                            ATC=AT*Czz{r00};
                        else
                            ATC=A_times_k_kron_B(AT,hz{r00},oo)+AT*Eu{oo};
                        end
                    elseif oo==3
                        matsizes=siz.nz*ones(3,2);
                        perms3={[2,1,3],[2,3,1]};
                        if precompute_czzzzz
                            if isempty(Czzz{r00})
                                Czzz{r00}=Pfunc(kron(hz{r00},Eu{2}),...
                                    matsizes,perms3{:})+...
                                    kronall(hz{r00},hz{r00},hz{r00});
                            end
                            ATC=AT*Czzz{r00};
                        else
                            ATC=A_times_k_kron_B(AT,hz{r00},oo);
                            ATC=ATC+A_times_sum_perms(AT,{hz{r00},Eu{2}},...
                                matsizes,true,perms3{:});
                        end
                    elseif oo==4
                        matsizes=siz.nz*ones(4,2);
                        % the following seems faster on smaller models
                        % and also potentially more accurate. Since we are
                        % to compute it several times, let's store it to
                        % save time instead of doing the alternative
                        %--------------------------------------------------
                        perms4={[1,3,2,4],[1,3,4,2],[3,1,4,2],[3,4,1,2],...
                            [3,1,2,4]};
                        if precompute_czzzzz
                            if isempty(Czzzz{r00})
                                Czzzz{r00}=kronall(hz{r00},hz{r00},hz{r00},hz{r00})+...
                                    Pfunc(kronall(hz{r00},hz{r00},Eu{2}),...
                                    matsizes,perms4{:})+Eu{oo};
                            end
                            ATC=AT*Czzzz{r00};
                        else
                            ATC=A_times_k_kron_B(AT,hz{r00},oo)+AT*Eu{oo};
                            ATC=ATC+A_times_sum_perms(AT,{hz{r00},hz{r00},Eu{2}},...
                                matsizes,true,perms4{:});
                        end
                    elseif oo==5
                        matsizes=siz.nz*ones(5,2);
                        perms51={[1,2,4,3,5],[1,2,4,5,3],[1,4,2,5,3],...
                            [1,4,5,2,3],[4,1,5,2,3],[4,5,1,2,3],...
                            [4,1,2,5,3],[4,1,2,3,5],[1,4,2,3,5]};
                        perms52={[2,1,3,4,5],[2,3,1,4,5],[2,3,4,1,5],...
                            [2,3,4,5,1]};
                        if precompute_czzzzz
                            if isempty(Czzzzz{r00})
                                Czzzzz{r00}=Pfunc(kronall(hz{r00},hz{r00},hz{r00},Eu{2}),...
                                    matsizes,perms51{:})+...
                                    Pfunc(kron(hz{r00},Eu{4}),matsizes,perms52{:})+...
                                    kronall(hz{r00},hz{r00},hz{r00},hz{r00},hz{r00});
                            end
                            ATC=AT*Czzzzz{r00};
                        else
                            ATC=A_times_k_kron_B(AT,hz{r00},oo);
                            ATC=ATC+A_times_sum_perms(AT,{hz{r00},hz{r00},hz{r00},Eu{2}},...
                                matsizes,true,perms51{:});
                            ATC=ATC+A_times_sum_perms(AT,{hz{r00},Eu{4}},...
                                matsizes,true,perms52{:});
                        end
                    elseif oo>5
                        error('perturbations implemented only up to order 5')
                    end
                    if accelerate
                        ATC=ATC(:,shrink{oo});
                    end
                end
            end
        end
    end
end


function Y=fvvv_vx_vx_vx(dvvv,vz)
[nd]=size(dvvv,1);
[~,nz]=size(vz);

if has_non_zeros(dvvv,vz)
    Y=utils.kronecker.A_times_k_kron_B(dvvv,vz,3);
else
    Y=zeros(nd,nz^3);
end

end

function flag=has_non_zeros(varargin)
flag=utils.cr.is_computable(varargin{:});
end

function Y=fvv_vx_vxx_omega_1(dvv,vz,vzz)
[~,nz]=size(vz);
if has_non_zeros(dvv,vz,vzz)
    Y=utils.kronecker.A_times_B_kron_C(dvv,vz,vzz);
    Y = utils.cr.dv_vz_omega(Y,nz,1);
else
    nd=size(dvv,1);
    Y=zeros(nd,nz^3);
end
end

