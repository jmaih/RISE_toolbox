function [y,newp,retcode]=sstate_model(obj,y,p,d,id) %#ok<INUSL>
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

retcode=0;
if nargin==1
    
    y={'LG_C','LG_I','LG_H','Z_C','Z_I','ZL_C','ZL_I','ZG_C','ZG_I','AL',...
        'AG','A','G_C','G_I','G_H','H_C','K_C','K_I','I_C','I_I','I','H_I',...
        'H','C','LAMBDA_C','LAMBDA_I','XI_C','XI_I'};

	% list of parameters herein calculated
	%-------------------------------------
	newp={};
	
    % flags on the calculation
    %--------------------------
    retcode=struct('unique',true,'imposed',true);
else
    % parameters to push
    %------------------
    newp=[];
    
    newp.theta_I=p.theta_C;
    
	AL=1;

	AG=d.ag;
   
	A=AL*AG;

    G_C = AG*(d.zg_I^p.theta_C)*(d.zg_C^(1-p.theta_C));
    
    G_I = AG*d.zg_I;
    
    G_H = AG;
    
    H_C = (1-p.theta_C)*((G_C-p.beta*p.gam)/(G_C-p.gam))*AG;
    
    kcssa = (1-newp.theta_I)*AG*d.zg_I;
    
    kcssb = (G_C-p.beta*p.gam)/(G_C-p.gam);
    
    kcssc = p.theta_C*(AG*d.zg_I-p.beta*(1-p.delta_I));
    
    kcssd = newp.theta_I*(AG*d.zg_I-p.beta*(1-p.delta_C));
    
    kcsse = ((p.beta*newp.theta_I)/(AG*d.zg_I-p.beta*(1-p.delta_I)))^(1/(1-newp.theta_I));
    
    K_C = kcssa*kcssb*(kcssc/kcssd)*kcsse;
    
    kissa = p.beta*newp.theta_I*(AG*d.zg_I-1+p.delta_C);
    
    kissb = AG*d.zg_I - p.beta*(1-p.delta_I);
    
    kissc = p.beta*newp.theta_I*(AG*d.zg_I-1+p.delta_I);
    
    K_I = (kissa/(kissb-kissc))*K_C;
    
    I_C = (AG*d.zg_I-1+p.delta_C)*K_C;
    
    I_I = (AG*d.zg_I-1+p.delta_I)*K_I;
    
    I = I_C + I_I;
    
    H_I = (1/d.zg_I)*((I/(K_I^newp.theta_I))^(1/(1-newp.theta_I)));
    
    H = H_C + H_I;
    
    C = (K_C^p.theta_C)*((d.zg_C*H_C)^(1-p.theta_C));
    
    LAMBDA_C = ((G_C-p.beta*p.gam)/(G_C-p.gam))/C;
    
    LAMBDA_I = H_I/((1-newp.theta_I)*AG*I);
    
    XI_C = LAMBDA_I;
    
    XI_I = XI_C;
    
	ZL_C=1;

	ZL_I=1;

	ZG_C=d.zg_C;
    
	ZG_I=d.zg_I;
    
    Z_C=ZL_C*ZG_C;
	
    Z_I=ZL_I*ZG_I;
    
    LG_C=log(G_C);
    
    LG_I=log(G_I);
    
    LG_H=log(G_H);
	
    newp.eta_C=AG;
    
    newp.eta_I=AG; 
    
    newp.kappa_C=I_C/K_C; 
    
    newp.kappa_I=I_I/K_I;
    
    ys=[LG_C,LG_I,LG_H,Z_C,Z_I,ZL_C,ZL_I,ZG_C,ZG_I,AL,AG,A,G_C,G_I,G_H,H_C,...
        K_C,K_I,I_C,I_I,I,H_I,H,C,LAMBDA_C,LAMBDA_I,XI_C,XI_I];
    
    ys=ys(:);
    
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

end