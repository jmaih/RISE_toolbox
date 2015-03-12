function [ys,obj,retcode,imposed]=sstate_model(obj,flag)

persistent yss

retcode=0;
imposed=false;
if flag==0
    ys={'X_H','X_F','TAU_H','TAU_F','P_H','D_H','D_F','R','Q_H','Q_F',...
        'P_F','P_A','P_B','Y_A','Y_B','K_H','K_F','L_H','L_F','W_H','W_F',...
        'I_H','I_F','C_H','C_F','A_H','B_H','A_F','B_F','LAMBDA_H',...
        'LAMBDA_F','XI_H','XI_F','N_H','N_F','RER','TOT','CTILDE_H',...
        'CTILDE_F','ITILDE_H','ITILDE_F','G_CH','G_IH','R_GC_H','R_CFH',...
        'R_IFH','R_GC_F',...
        'G_F','G_H','M_F','M_H','V_F','V_H','V_HF','Z_F','Z_H','Z_HF'};
    % initialize the persistent variable at first call
    %--------------------------------------------------
    yss=[];
else
    if isempty(yss)
        p=get(obj,'parameters');
        newp=struct();
        if p.theta == 1
            p.theta = 1.000001;
            newp.theta=p.theta;
        elseif p.theta < 0.01
            p.theta  = 0.01;
            newp.theta=p.theta;
        end
        
        v=p.vss_H;
        z=p.zss_H;
        
        G_F=p.gss_F;
        G_H=p.gss_H;
        M_F=p.mss_F;
        M_H=p.mss_H;
        V_F=v;
        V_H=v;
        V_HF=p.v_hf_ss;
        Z_F=z;
        Z_H=z;
        Z_HF=p.z_hf_ss;
        
        X_H = 1/v;
        X_F = 1/v;
        
        TAU_H = p.gss_H;
        TAU_F = p.gss_F;
        
        P_H = 1;
        
        D_H = 0;
        D_F = 0;
        
        R = (v^(p.alpha/(1-p.alpha))*z)^(1-p.mu*(1-p.gam))/p.beta;
        
        Q_H =  R - (1-p.delta)/v;
        Q_F =  R - (1-p.delta)/v;
        
        P_F = fsolve(@pfssfn,1,optimset('Display','none'));
        
        [P_A,P_B,Y_A,Y_B,A_H,A_F]=many_variables(P_F);
        
        K_H = (p.alpha/Q_H)*P_A*Y_A;
        K_F = (p.alpha/Q_F)*(P_B/P_F)*Y_B;
        
        L_H = (1/z)*((Q_H/p.alpha)*(1/P_A))^(p.alpha/(1-p.alpha))*Y_A;
        L_F = (1/z)*((Q_F/p.alpha)*(P_F/P_B))^(p.alpha/(1-p.alpha))*Y_B;
        
        W_H = (1-p.alpha)*z*(p.alpha/Q_H)^(p.alpha/(1-p.alpha))*P_A^(1/(1-p.alpha));
        W_F = (1-p.alpha)*z*(p.alpha/Q_F)^(p.alpha/(1-p.alpha))*(P_B/P_F)^(1/(1-p.alpha));
        
        newp.eta_H = v^(1/(1-p.alpha))*z - 1 + p.delta;
        newp.eta_F = v^(1/(1-p.alpha))*z - 1 + p.delta;
        
        I_H = newp.eta_H*K_H;
        I_F = newp.eta_F*K_F;
        
        C_H = (1-p.alpha)*(p.mu/(1-p.mu))*(z*p.mss_H*(p.alpha/Q_H)^(p.alpha/(1-p.alpha))*P_A^(1/(1-p.alpha))-P_A*Y_A);
        C_F = (1-p.alpha)*(p.mu/(1-p.mu))*(z*p.mss_F*(p.alpha/Q_F)^(p.alpha/(1-p.alpha))*(P_B/P_F)^(1/(1-p.alpha))-(P_B/P_F)*Y_B);
        
        B_H = p.omega*Y_A*P_A*P_B^(-p.theta);
        
        B_F = (1-p.omega)*Y_B*(P_B/P_F)^(1-p.theta);
        
        LAMBDA_H = p.mu*((C_H^p.mu*(1-L_H/p.mss_H)^(1-p.mu))^(1-p.gam))/C_H;
        LAMBDA_F = p.mu*((C_F^p.mu*(1-L_F/p.mss_F)^(1-p.mu))^(1-p.gam))/C_F;
        
        XI_H = LAMBDA_H/v;
        XI_F = LAMBDA_F/v;
        
        N_H = P_A*A_F/(p.v_hf_ss^(p.alpha/(1-p.alpha))*p.z_hf_ss)-P_B*B_H;
        N_F = (P_B*B_H*p.v_hf_ss^(p.alpha/(1-p.alpha))*p.z_hf_ss-P_A*A_F)/P_F;
        
        RER = P_F;
        
        TOT = P_B/P_A;
        
        CTILDE_H = C_H + p.gss_H;
        CTILDE_F = C_F + p.gss_F;
        
        ITILDE_H = I_H;
        ITILDE_F = I_F;
        
        G_CH = v^(p.alpha/(1-p.alpha))*z;
        G_IH = v^(1/(1-p.alpha))*z;
        R_GC_H = p.gss_H/C_H;
        R_CFH = (C_F/C_H)/(p.v_hf_ss^(p.alpha/(1-p.alpha))*p.z_hf_ss);
        R_IFH = (I_F/I_H)/(p.v_hf_ss^(1/(1-p.alpha))*p.z_hf_ss);
        R_GC_F = p.gss_F/C_F;
        
        ys=[X_H,X_F,TAU_H,TAU_F,P_H,D_H,D_F,R,Q_H,Q_F,P_F,P_A,P_B,Y_A,Y_B,K_H,...
            K_F,L_H,L_F,W_H,W_F,I_H,I_F,C_H,C_F,A_H,B_H,A_F,B_F,LAMBDA_H,...
            LAMBDA_F,XI_H,XI_F,N_H,N_F,RER,TOT,CTILDE_H,CTILDE_F,ITILDE_H,...
            ITILDE_F,G_CH,G_IH,R_GC_H,R_CFH,R_IFH,R_GC_F,...
            G_F,G_H,M_F,M_H,V_F,V_H,V_HF,Z_F,Z_H,Z_HF].';
        
        % push the new parameters
        %-------------------------
        obj=set(obj,'parameters',newp);
        
        % populate the persistent variable
        %---------------------------------
        yss=ys;
    else
        % do not recompute the series. Just load the unchanged past values
        %------------------------------------------------------------------
        ys=yss;
    end
end

    function resid=pfssfn(pf)
        [~,~,Y_A,~,A_H,A_F]=many_variables(pf);
        resid = Y_A - A_H - A_F/(p.v_hf_ss^(p.alpha/(1-p.alpha))*p.z_hf_ss);
    end

    function [P_A,P_B,Y_A,Y_B,A_H,A_F]=many_variables(P_F)
        P_A = ((1-p.omega-p.omega*P_F^(1-p.theta))/(1-2*p.omega))^(1/(1-p.theta));
        P_B = (((1-p.omega)*P_F^(1-p.theta)-p.omega)/(1-2*p.omega))^(1/(1-p.theta));
        Y_A = (p.mu*(1-p.alpha)*z*p.mss_H*(p.alpha*P_A/Q_H)^(p.alpha/(1-p.alpha))+(1-p.mu)*p.gss_H/P_A) ...
            /(1-p.alpha*p.mu-(1/v)*(1-p.mu)*(v^(1/(1-p.alpha))*z-1+p.delta)*p.alpha/Q_H);
        
        Y_B = (p.mu*(1-p.alpha)*z*p.mss_F*((p.alpha/Q_F)*(P_B/P_F))^(p.alpha/(1-p.alpha))+(1-p.mu)*P_F*p.gss_F/P_B) ...
            /(1-p.alpha*p.mu-(1/v)*(1-p.mu)*(v^(1/(1-p.alpha))*z-1+p.delta)*p.alpha/Q_F);
        A_H = (1-p.omega)*Y_A*P_A^(1-p.theta);
        A_F = p.omega*Y_B*(P_B/P_F)*(P_A/P_F)^(-p.theta);
    end
end