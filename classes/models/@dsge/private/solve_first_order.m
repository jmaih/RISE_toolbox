
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
            df_plus=dv{rt,rplus}(:,pos.v.f_plus);
            db_plus=dv{rt,rplus}(:,pos.v.b_plus);
            df_Lf_Tzp_Lp(:,pos.t.p)=df_plus*Tz_pb(pos.t.f,1:siz.np,rplus);% place in the p position
            db_Lb_Tzb_Lb(:,pos.t.b)=db_plus*Tz_pb(pos.t.b,siz.np+(1:siz.nb),rplus);% place in the b position
            A0sig(:,:,rt) = A0sig(:,:,rt) + A0_0_1 + df_Lf_Tzp_Lp + db_Lb_Tzb_Lb;
            
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
    % ensure the result is sparse
    %-----------------------------
    for rt=1:siz.h
        Tz{rt}=sparse(Tz{rt});
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