function [y,newp,retcode]=sstate_model12(obj,y,p,d,id) %#ok<INUSL>
% sstate_model -- shows the way of writing a RISE steady state file
%
% ::
%
%
%   [y,newp,retcode]=sstate_model(obj,y,p,d,id)
%
% Args:
%
%    - **obj** [rise|dsge]: model object (not always needed)
%
%    - **y** [vector]: endo_nbr x 1 vector of initial steady state
%
%    - **p** [struct]: parameter structure
%
%    - **d** [struct]: definitions
%
%    - **id** [vector]: location of the variables to calculate
%
% Returns:
%    :
%
%      CASE 1: one input argument
%
%    - **y** [cell array]: list of the variables for which the steady state
%    will be calculated within the steady state function
%
%    - **newp** [cell array]: List of the parameters calculated inside the
%    steady state function
%
%    - **retcode** [struct]: possible fields are "imposed", "unique", "loop".
%    The default value for all of those is false.
%      - "imposed": This tells RISE not to check that this is actually solves
%          the steady state. Hence, RISE will attempt to approximate the model
%          around the chosen point
%      - "unique": This tells RISE that the steady state is the same across
%          all regimes. RISE will call the function only once but instead of
%          using just any parameter vector, it will use the ergodic mean of
%          the parameter vector (mean across all regimes).
%      - "loop": This tells RISE that if the user was given the steady state
%          values for some of the variables, he would be able to compute the
%          steady state for the remaining variables. RISE will then exploit
%          this information to optimize over the variables that the user needs
%          for computing the steady state.
%
%      CASE 2: More than one input argument
%
%    - **y** []: endo_nbr x 1 vector of updated steady state
%
%    - **newp** [struct]: structure containing updated parameters if any
%
%    - **retcode** [0|number]: return 0 if there are no problems, else return
%      any number different from 0
%
% Note:
%
%    - If the user knows the steady state, it is always an advantage. If the
%    steady state is computed numerically, we don't know whether it is unique
%    or not. Not that it really matters but... some economists have a strong
%    aversion towards models with multiple equilibria.
%
%    - If the user does not know the solution for all the variables in the
%    steady state, it is a good idea to take a log-linear approximation for
%    the variables that potentially have a nonzero steady state. Hence the
%    user should give that information to RISE.
%
%    - One can potentially improve on the above point by explicit bounds on
%    the variables. But this is not implemented.
%
%    - An alternative that potentially avoids taking a loglinearization is to
%    to reset the values proposed by the optimizer whenever they are in a bad
%    region. It is unclear whether this always works.
%
%    - So be on the safe side, i.e. don't do like me: compute your steady
%    state analytically.
%
% Example:
%
%    See also:

persistent yss newp_ ncp_shocks

retcode=0;
if nargin==1
    y={'X_H','X_F','TAU_H','TAU_F','P_H','D_H','D_F','R','Q_H','Q_F',...
        'P_F','P_A','P_B','Y_A','Y_B','K_H','K_F','L_H','L_F','W_H','W_F',...
        'I_H','I_F','C_H','C_F','A_H','B_H','A_F','B_F','LAMBDA_H',...
        'LAMBDA_F','XI_H','XI_F','N_H','N_F','RER','TOT','CTILDE_H',...
        'CTILDE_F','ITILDE_H','ITILDE_F','G_CH','G_IH','R_GC_H','R_CFH',...
        'R_IFH','R_GC_F',...
        'G_F','G_H','M_F','M_H','V_F','V_H','V_HF','Z_F','Z_H','Z_HF'};
		
    % flag the model of interest
    %---------------------------
    ncp_shocks=any(strcmp(obj.endogenous.name,'M_HF'));
    
    if ncp_shocks
        y=[y,{'G_LH','R_LFH','M_HF'}];
    end
	
    % initialize the persistent variable at first call
    %--------------------------------------------------
    yss=[];

	% list of parameters herein calculated
	%--------------------------------------
	newp={'theta','eta_H','eta_F'};
	
    % flags on the calculation
    %--------------------------
    retcode=struct('unique',true,'imposed',true);
else
    if isempty(yss)
        
        newp=struct();
        
        if p.theta == 1
            
            p.theta = 1.000001;
            
        elseif p.theta < 0.01
            
            p.theta  = 0.01;
            
        end
        
        newp.theta=p.theta;
        
        v=p.vss_H;
        z=p.zss_H;
        
        G_F=p.gss_F;
        G_H=p.gss_H;
        m=p.mss_H;
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
        
        if ncp_shocks
            M_HF=p.m_hf_ss;
            R = (m*v^(p.alpha/(1-p.alpha))*z)^(1-p.mu*(1-p.gam))/p.beta;
        else
            R = (v^(p.alpha/(1-p.alpha))*z)^(1-p.mu*(1-p.gam))/p.beta;
        end
        
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
        
        if ncp_shocks
            newp.eta_H = m*v^(1/(1-p.alpha))*z - 1 + p.delta;
            newp.eta_F = m*v^(1/(1-p.alpha))*z - 1 + p.delta;
        else
            newp.eta_H = v^(1/(1-p.alpha))*z - 1 + p.delta;
            newp.eta_F = v^(1/(1-p.alpha))*z - 1 + p.delta;
        end
        
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
        
        if ncp_shocks
            N_H = P_A*A_F/(M_HF*p.v_hf_ss^(p.alpha/(1-p.alpha))*p.z_hf_ss)-P_B*B_H;
            N_F = (P_B*B_H*M_HF*p.v_hf_ss^(p.alpha/(1-p.alpha))*p.z_hf_ss-P_A*A_F)/P_F;
        else
            N_H = P_A*A_F/(p.v_hf_ss^(p.alpha/(1-p.alpha))*p.z_hf_ss)-P_B*B_H;
            N_F = (P_B*B_H*p.v_hf_ss^(p.alpha/(1-p.alpha))*p.z_hf_ss-P_A*A_F)/P_F;
        end
        
        RER = P_F;
        
        TOT = P_B/P_A;
        
        CTILDE_H = C_H + p.gss_H;
        CTILDE_F = C_F + p.gss_F;
        
        ITILDE_H = I_H;
        ITILDE_F = I_F;
        
        G_CH = v^(p.alpha/(1-p.alpha))*z;
        G_IH = v^(1/(1-p.alpha))*z;
        if ncp_shocks
            G_CH=m*G_CH;
            G_IH=m*G_IH;
            G_LH=m;
            R_LFH=L_F/L_H*1/M_HF;
        end
        R_GC_H = p.gss_H/C_H;
        R_CFH = (C_F/C_H)/(M_HF*p.v_hf_ss^(p.alpha/(1-p.alpha))*p.z_hf_ss);
        R_IFH = (I_F/I_H)/(M_HF*p.v_hf_ss^(1/(1-p.alpha))*p.z_hf_ss);
        R_GC_F = p.gss_F/C_F;
        
        ys=[X_H,X_F,TAU_H,TAU_F,P_H,D_H,D_F,R,Q_H,Q_F,P_F,P_A,P_B,Y_A,Y_B,K_H,...
            K_F,L_H,L_F,W_H,W_F,I_H,I_F,C_H,C_F,A_H,B_H,A_F,B_F,LAMBDA_H,...
            LAMBDA_F,XI_H,XI_F,N_H,N_F,RER,TOT,CTILDE_H,CTILDE_F,ITILDE_H,...
            ITILDE_F,G_CH,G_IH,R_GC_H,R_CFH,R_IFH,R_GC_F,...
            G_F,G_H,M_F,M_H,V_F,V_H,V_HF,Z_F,Z_H,Z_HF];
        if ncp_shocks
            ys=[ys,G_LH,R_LFH,M_HF];
        end
        
        ys=ys(:);
        
        % populate the persistent variable
        %---------------------------------
        yss=ys;
        
        newp_=newp;
        
    else
        % do not recompute the series. Just load the unchanged past values
        %------------------------------------------------------------------
        ys=yss;
        
        newp=newp_;
        
    end
    % check the validity of the calculations
    %----------------------------------------
    if ~utils.error.valid(ys)
        
        retcode=1;
        
    else
        % push the calculations
        %----------------------
        y(id)=ys;
        
    end
    
end

    function resid=pfssfn(pf)
        [~,~,Y_A,~,A_H,A_F]=many_variables(pf);
        resid = Y_A - A_H - A_F/(p.v_hf_ss^(p.alpha/(1-p.alpha))*p.z_hf_ss);
    end

    function [P_A,P_B,Y_A,Y_B,A_H,A_F]=many_variables(P_F)
        P_A = ((1-p.omega-p.omega*P_F^(1-p.theta))/(1-2*p.omega))^(1/(1-p.theta));
        P_B = (((1-p.omega)*P_F^(1-p.theta)-p.omega)/(1-2*p.omega))^(1/(1-p.theta));
        if ncp_shocks
            Y_A = (p.mu*(1-p.alpha)*z*p.mss_H*(p.alpha*P_A/Q_H)^(p.alpha/(1-p.alpha))+(1-p.mu)*p.gss_H/P_A) ...
                /(1-p.alpha*p.mu-(1/v)*(1-p.mu)*(m*v^(1/(1-p.alpha))*z-1+p.delta)*p.alpha/Q_H);
            
            Y_B = (p.mu*(1-p.alpha)*z*p.mss_F*((p.alpha/Q_F)*(P_B/P_F))^(p.alpha/(1-p.alpha))+(1-p.mu)*P_F*p.gss_F/P_B) ...
                /(1-p.alpha*p.mu-(1/v)*(1-p.mu)*(m*v^(1/(1-p.alpha))*z-1+p.delta)*p.alpha/Q_F);
        else
            Y_A = (p.mu*(1-p.alpha)*z*p.mss_H*(p.alpha*P_A/Q_H)^(p.alpha/(1-p.alpha))+(1-p.mu)*p.gss_H/P_A) ...
                /(1-p.alpha*p.mu-(1/v)*(1-p.mu)*(v^(1/(1-p.alpha))*z-1+p.delta)*p.alpha/Q_H);
            
            Y_B = (p.mu*(1-p.alpha)*z*p.mss_F*((p.alpha/Q_F)*(P_B/P_F))^(p.alpha/(1-p.alpha))+(1-p.mu)*P_F*p.gss_F/P_B) ...
                /(1-p.alpha*p.mu-(1/v)*(1-p.mu)*(v^(1/(1-p.alpha))*z-1+p.delta)*p.alpha/Q_F);
        end
        A_H = (1-p.omega)*Y_A*P_A^(1-p.theta);
        A_F = p.omega*Y_B*(P_B/P_F)*(P_A/P_F)^(-p.theta);
    end

end