function [obj,retcode]=second_order_solver(obj,syst_mat)

Defaults=struct();
if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    obj=Defaults;
    return
end
endo_nbr=obj.endogenous.number(end);
exo_nbr=sum(obj.exogenous.number);
h=obj.markov_chains.regimes_number;
shock_horizon=max(obj.exogenous.shock_horizon);
if shock_horizon>1 
    error('RISE does not currently solve up to a 2nd order models in which agents anticipate the future')
end

if obj.is_hybrid_expectations_model||obj.is_sticky_information_model||...
        obj.is_optimal_policy_model
    error('RISE does not currently solve optimal policy, hybrid expectations or sticky information models up to a 2nd order')
end

theta_hat=obj.solution.theta_hat;
solve_disable_theta=obj.options.solve_disable_theta;
for ireg=1:h
    if isempty(theta_hat{ireg})
        break
    end
    theta_hat{ireg}=(1-solve_disable_theta)*theta_hat{ireg};
end
[Gplus01,A0]=msre_solvers.integrate_structure(syst_mat,h);

A0Gm_inv=msre_solvers.compute_preconditioner(Gplus01,A0,obj.solution.Q,obj.solution.m_x);
GBAR=cell(h);
for s0=1:h
    for s1=1:h
        GBAR{s0,s1}=sparse(A0Gm_inv{s0}*Gplus01{s0,s1});
    end
end
clear Gplus01

Ie=speye(exo_nbr);
% Ix=speye(endo_nbr);
optim_opt=obj.options;
if  obj.options.debug
    optim_opt.fix_point_verbose=true;
end

% solving for m_x_x
%------------------
[Dxx_bar,m_xx_func]=routine_mxx();

m_x_x=zeros(endo_nbr^3*h,1);
[m_x_x,retcode,tau]=transpose_free_quasi_minimum_residual(...
    m_xx_func,... % coefficient matrix
    -Dxx_bar(:),... % right hand side
    m_x_x,... % initial guess
    optim_opt.fix_point_TolFun,... % tolerance level
    optim_opt.fix_point_maxiter,... % maximum number of iterations
    optim_opt.fix_point_verbose); % flag for printing progress or not
if retcode==201 && tau<=optim_opt.fix_point_TolFun
    retcode=0;
end
if retcode
    return
end
m_x_x=reshape(m_x_x,endo_nbr,endo_nbr^2,h);
clear Dxx
% solving for m_x_e
%------------------
obj.solution.m_x_e=routine_mxe();

% solving for m_x_sig
%--------------------
obj.solution.m_x_sig=routine_mxsig();

% solving for m_e_e
%--------------------
obj.solution.m_e_e=routine_mee();

% solving for m_e_sig
%--------------------
obj.solution.m_e_sig=routine_mesig();

% solving for m_sig_sig
%----------------------
[Dsigsig_bar,m_sig_sig_func]=routine_msigsig();

m_sig_sig=zeros(endo_nbr*1*h,1);
[m_sig_sig,retcode,tau]=transpose_free_quasi_minimum_residual(...
    m_sig_sig_func,... % coefficient matrix
    -Dsigsig_bar(:),... % right hand side
    m_sig_sig,... % initial guess
    optim_opt.fix_point_TolFun,... % tolerance level
    optim_opt.fix_point_maxiter,... % maximum number of iterations
    optim_opt.fix_point_verbose); % flag for printing progress or not
if retcode==201 && tau<=optim_opt.fix_point_TolFun
    retcode=0;
end
if retcode
    return
end
m_sig_sig=reshape(m_sig_sig,endo_nbr,1,h);

% store further output
%---------------------
for s0=1:h
    obj.solution.m_x_x{s0}=sparse(m_x_x(:,:,s0));
    obj.solution.m_sig_sig{s0}=sparse(m_sig_sig(:,:,s0));
end
clear m_x_x m_sig_sig

    function [Dsigsig_bar,AB]=routine_msigsig()
        Dsigsig_bar=zeros(endo_nbr,1,h);
        for s00=1:h
            Dsigsig=0;
            ms0ms0=kron(obj.solution.m_sig{s00},obj.solution.m_sig{s00});
            for s11=1:h
                mx1ms0=obj.solution.m_x{s11}*obj.solution.m_sig{s00};
                mxs1ms0=obj.solution.m_x_sig{s11}*obj.solution.m_sig{s00};
                Dsigsig=Dsigsig+...
                    syst_mat.Gpp{s00,s11}*kron(mx1ms0+obj.solution.m_sig{s11},mx1ms0)+...
                    syst_mat.Gpc{s00,s11}*kron(obj.solution.m_sig{s00},mx1ms0)+...
                    syst_mat.Gp{s00,s11}*(m_x_x(:,:,s11)*ms0ms0+mxs1ms0)+...
                    (syst_mat.Gpp{s00,s11}*kron(obj.solution.m_e{s11},obj.solution.m_e{s11})+...
                    syst_mat.Gp{s00,s11}*obj.solution.m_e_e{s11})*Ie(:)+...
                    syst_mat.Gpp{s00,s11}*kron(mx1ms0+obj.solution.m_sig{s11},obj.solution.m_sig{s11})+...
                    syst_mat.Gpc{s00,s11}*kron(obj.solution.m_sig{s00},obj.solution.m_sig{s11})+...
                    syst_mat.Gcp{s00,s11}*kron(mx1ms0,obj.solution.m_sig{s00})+...
                    syst_mat.Gcc{s00,s11}*ms0ms0+...
                    syst_mat.Gtt{s00,s11}*kron(obj.solution.theta_hat{s11},obj.solution.theta_hat{s11})+...
                    syst_mat.Gp{s00,s11}*mxs1ms0;
            end
            Dsigsig_bar(:,:,s00)=-sparse(A0Gm_inv{s00}*Dsigsig);
        end
        AB=@matrix_vector_product;
        function lhs=matrix_vector_product(m_sig_sig)
            m_sig_sig=reshape(m_sig_sig,endo_nbr,1,h);
            lhs=zeros(endo_nbr,1,h);
            for s000=1:h
                for s111=1:h
                    lhs(:,:,s000)=lhs(:,:,s000)+GBAR{s000,s111}*m_sig_sig(:,:,s111);
                end
                lhs(:,:,s000)=lhs(:,:,s000)+m_sig_sig(:,:,s000);
            end
            lhs=lhs(:);
        end
    end

    function mesig=routine_mesig()
        mesig=cell(1,h);
        for s00=1:h
            Desig=0;
            mse00=kron(obj.solution.m_sig{s00},obj.solution.m_e{s00});
            for s11=1:h
                mxmsms=obj.solution.m_x{s11}*obj.solution.m_sig{s00}+obj.solution.m_sig{s11};
                mx1me0=obj.solution.m_x{s11}*obj.solution.m_e{s00};
                Desig=Desig+...
                    syst_mat.Gpp{s00,s11}*kron(mxmsms,mx1me0)+...
                    syst_mat.Gpc{s00,s11}*kron(obj.solution.m_sig{s00},mx1me0)+...
                    syst_mat.Gpt{s00,s11}*kron(obj.solution.theta_hat{s11},mx1me0)+...
                    syst_mat.Gp{s00,s11}*obj.solution.m_x_sig{s11}*obj.solution.m_e{s00}+...
                    syst_mat.Gcp{s00,s11}*kron(mxmsms,obj.solution.m_e{s00})+...
                    syst_mat.Gcc{s00,s11}*mse00+...
                    syst_mat.Gct{s00,s11}*kron(obj.solution.theta_hat{s11},obj.solution.m_e{s00})+...
                    utils.kronecker.A_times_kron_B_I(syst_mat.Gep{s00,s11},mxmsms,exo_nbr)+...syst_mat.Gep{s00,s11}*kron(mxmsms,Ie)+...
                    utils.kronecker.A_times_kron_B_I(syst_mat.Gec{s00,s11},obj.solution.m_sig{s00},exo_nbr)+...syst_mat.Gec{s00,s11}*kron(obj.solution.m_sig{s00},Ie)+...
                    utils.kronecker.A_times_kron_B_I(syst_mat.Get{s00,s11},obj.solution.theta_hat{s11},exo_nbr);%syst_mat.Get{s00,s11}*kron(obj.solution.theta_hat{s11},Ie);
            end
            mesig{s00}=-sparse(A0Gm_inv{s00}*Desig);
        end
    end

    function mee=routine_mee()
        mee=cell(1,h);
        for s00=1:h
            Dee=0;
            for s11=1:h
                mxe10=obj.solution.m_x{s11}*obj.solution.m_e{s00};
                Dee=Dee+...
                    syst_mat.Gpp{s00,s11}*kron(mxe10,mxe10)+...
                    syst_mat.Gpc{s00,s11}*kron(obj.solution.m_e{s00},mxe10)+...
                    utils.kronecker.A_times_kron_I_B(syst_mat.Gpe{s00,s11},mxe10,exo_nbr)+...syst_mat.Gpe{s00,s11}*kron(Ie,mxe10)+...
                    syst_mat.Gp{s00,s11}*m_x_x(:,:,s11)*kron(obj.solution.m_e{s00},obj.solution.m_e{s00})+...
                    syst_mat.Gcp{s00,s11}*kron(mxe10,obj.solution.m_e{s00})+...
                    syst_mat.Gcc{s00,s11}*kron(obj.solution.m_e{s00},obj.solution.m_e{s00})+...
                    utils.kronecker.A_times_kron_I_B(syst_mat.Gce{s00,s11},obj.solution.m_e{s00},exo_nbr)+...syst_mat.Gce{s00,s11}*kron(Ie,obj.solution.m_e{s00})+...
                    utils.kronecker.A_times_kron_B_I(syst_mat.Gep{s00,s11},mxe10,exo_nbr)+...syst_mat.Gep{s00,s11}*kron(mxe10,Ie)+...
                    utils.kronecker.A_times_kron_B_I(syst_mat.Gec{s00,s11},obj.solution.m_e{s00},exo_nbr)+...syst_mat.Gec{s00,s11}*kron(obj.solution.m_e{s00},Ie)+...
                    syst_mat.Gee{s00,s11};
            end
            mee{s00}=-sparse(A0Gm_inv{s00}*Dee);
        end
    end

    function mxsig=routine_mxsig()
        mxsig=cell(1,h);
        for s00=1:h
            Dxsig=0;
            msigx00=kron(obj.solution.m_sig{s00},obj.solution.m_x{s00});
            for s11=1:h
                mx1s0=obj.solution.m_x{s11}*obj.solution.m_sig{s00};
                mx1x0=obj.solution.m_x{s11}*obj.solution.m_x{s00};
                Dxsig=Dxsig+...
                    syst_mat.Gpp{s00,s11}*kron(mx1s0,mx1x0)+...
                    syst_mat.Gpp{s00,s11}*kron(obj.solution.m_sig{s11},mx1x0)+...
                    syst_mat.Gpc{s00,s11}*kron(obj.solution.m_sig{s00},mx1x0)+...
                    syst_mat.Gpt{s00,s11}*kron(theta_hat{s11},mx1x0)+...
                    syst_mat.Gp{s00,s11}*m_x_x(:,:,s11)*msigx00+...
                    syst_mat.Gcp{s00,s11}*kron(mx1s0,obj.solution.m_x{s00})+...
                    syst_mat.Gcp{s00,s11}*kron(obj.solution.m_sig{s11},obj.solution.m_x{s00})+...
                    syst_mat.Gcc{s00,s11}*msigx00+...
                    utils.kronecker.A_times_kron_B_I(syst_mat.Gmp{s00,s11},mx1s0,endo_nbr)+...syst_mat.Gmp{s00,s11}*kron(mx1s0,Ix)+...
                    utils.kronecker.A_times_kron_B_I(syst_mat.Gmp{s00,s11},obj.solution.m_sig{s11},endo_nbr)+...syst_mat.Gmp{s00,s11}*kron(obj.solution.m_sig{s11},Ix)+...
                    utils.kronecker.A_times_kron_B_I(syst_mat.Gmc{s00,s11},obj.solution.m_sig{s00},endo_nbr);%syst_mat.Gmc{s00,s11}*kron(obj.solution.m_sig{s00},Ix);
            end
            mxsig{s00}=-sparse(A0Gm_inv{s00}*Dxsig);
        end
    end

    function mxe=routine_mxe()
        % N.B: There is something potentially confusing about the
        % computations done here. These results are obtained from
        % differentiating G with respect to e and then with respect to x
        % and not the other way around as one would have expected. This is
        % done so that we remain consistent with the kronecker products
        % xe=[x1_e1,x1_e2,...,x1_ek,x2_e1,x2_e2,...,x2_ek,...,
        % xn_e1,xn_e2,...,xn_ek]. 
        mxe=cell(1,h);
        for s00=1:h
            Dxe=0;
            for s11=1:h
                mx1x0=obj.solution.m_x{s11}*obj.solution.m_x{s00};
                mx1e0=obj.solution.m_x{s11}*obj.solution.m_e{s00};
                Dxe=Dxe+...
                    syst_mat.Gpp{s00,s11}*kron(mx1x0,mx1e0)+...
                    syst_mat.Gpc{s00,s11}*kron(obj.solution.m_x{s00},mx1e0)+...
                    utils.kronecker.A_times_kron_I_B(syst_mat.Gpm{s00,s11},mx1e0,endo_nbr)+...syst_mat.Gpm{s00,s11}*kron(Ix,mx1e0)+...
                    syst_mat.Gp{s00,s11}*m_x_x(:,:,s11)*kron(obj.solution.m_x{s00},obj.solution.m_e{s00})+...
                    syst_mat.Gcp{s00,s11}*kron(mx1x0,obj.solution.m_e{s00})+...
                    syst_mat.Gcc{s00,s11}*kron(obj.solution.m_x{s00},obj.solution.m_e{s00})+...
                    utils.kronecker.A_times_kron_I_B(syst_mat.Gcm{s00,s11},obj.solution.m_e{s00},endo_nbr)+...syst_mat.Gcm{s00,s11}*kron(Ix,obj.solution.m_e{s00})+...
                    utils.kronecker.A_times_kron_B_I(syst_mat.Gep{s00,s11},mx1x0,exo_nbr)+...syst_mat.Gep{s00,s11}*kron(mx1x0,Ie)+...
                    utils.kronecker.A_times_kron_B_I(syst_mat.Gec{s00,s11},obj.solution.m_x{s00},exo_nbr)+...syst_mat.Gec{s00,s11}*kron(obj.solution.m_x{s00},Ie)+...
                    syst_mat.Gem{s00,s11};
% % % %                 Dxe=Dxe+...
% % % %                     syst_mat.Gpp{s00,s11}*kron(obj.solution.m_x{s11}*obj.solution.m_e{s00},mx1x0)+...
% % % %                     syst_mat.Gpc{s00,s11}*kron(obj.solution.m_e{s00},mx1x0)+...
% % % %                     utils.kronecker.A_times_kron_I_B(syst_mat.Gpe{s00,s11},mx1x0,exo_nbr)+...syst_mat.Gpe{s00,s11}*kron(Ie,mx1x0)+...
% % % %                     syst_mat.Gp{s00,s11}*m_x_x(:,:,s11)*kron(obj.solution.m_e{s00},obj.solution.m_x{s00})+...
% % % %                     syst_mat.Gcp{s00,s11}*kron(obj.solution.m_x{s11}*obj.solution.m_e{s00},obj.solution.m_x{s00})+...
% % % %                     syst_mat.Gcc{s00,s11}*kron(obj.solution.m_e{s00},obj.solution.m_x{s00})+...
% % % %                     utils.kronecker.A_times_kron_I_B(syst_mat.Gce{s00,s11},obj.solution.m_x{s00},exo_nbr)+...syst_mat.Gce{s00,s11}*kron(Ie,obj.solution.m_x{s00})+...
% % % %                     utils.kronecker.A_times_kron_B_I(syst_mat.Gmp{s00,s11},obj.solution.m_x{s11}*obj.solution.m_e{s00},endo_nbr)+...syst_mat.Gmp{s00,s11}*kron(obj.solution.m_x{s11}*obj.solution.m_e{s00},Ix)+...
% % % %                     utils.kronecker.A_times_kron_B_I(syst_mat.Gmc{s00,s11},obj.solution.m_e{s00},endo_nbr)+...syst_mat.Gmc{s00,s11}*kron(obj.solution.m_e{s00},Ix)+...
% % % %                     syst_mat.Gme{s00,s11};
            end
            mxe{s00}=-sparse(A0Gm_inv{s00}*Dxe);
        end
    end

    function [Dxx_bar,AB]=routine_mxx()
        Dxx_bar=zeros(endo_nbr,endo_nbr^2,h);
        mx0x0=cell(1,h);
        m_x=obj.solution.m_x;
        for s00=1:h
            mx0x0{s00}=kron(m_x{s00},m_x{s00});
            for s11=1:h
                mx1x0=m_x{s11}*m_x{s00};
                Dxx_bar(:,:,s00)=Dxx_bar(:,:,s00)+utils.kronecker.A_times_k_kron_B(syst_mat.Gpp{s00,s11},mx1x0,2)+...syst_mat.Gpp{s00,s11}*kron(mx1x0,mx1x0)+...
                    syst_mat.Gpc{s00,s11}*kron(m_x{s00},mx1x0)+...
                    utils.kronecker.A_times_kron_I_B(syst_mat.Gpm{s00,s11},mx1x0,endo_nbr)+...syst_mat.Gpm{s00,s11}*kron(Ix,mx1x0)+... 
                    syst_mat.Gcp{s00,s11}*kron(mx1x0,m_x{s00})+...
                    syst_mat.Gcc{s00,s11}*mx0x0{s00}+...
                    utils.kronecker.A_times_kron_I_B(syst_mat.Gcm{s00,s11},m_x{s00},endo_nbr)+...syst_mat.Gcm{s00,s11}*kron(Ix,m_x{s00})+...
                    utils.kronecker.A_times_kron_B_I(syst_mat.Gmp{s00,s11},mx1x0,endo_nbr)+...syst_mat.Gmp{s00,s11}*kron(mx1x0,Ix)+...
                    utils.kronecker.A_times_kron_B_I(syst_mat.Gmc{s00,s11},m_x{s00},endo_nbr)+... syst_mat.Gmc{s00,s11}*kron(m_x{s00},Ix)+...
                    syst_mat.Gmm{s00,s11};
            end
            Dxx_bar(:,:,s00)=A0Gm_inv{s00}*Dxx_bar(:,:,s00);
        end
        AB=@matrix_vector_product;
        function lhs=matrix_vector_product(m_x_x)
            m_x_x=reshape(m_x_x,endo_nbr,endo_nbr^2,h);
            lhs=zeros(endo_nbr,endo_nbr^2,h);
            for s000=1:h
                for s111=1:h
                    lhs(:,:,s000)=lhs(:,:,s000)+GBAR{s000,s111}*m_x_x(:,:,s111);
                end
                lhs(:,:,s000)=lhs(:,:,s000)*mx0x0{s000}+m_x_x(:,:,s000);
            end
            lhs=lhs(:);
        end
    end
end

