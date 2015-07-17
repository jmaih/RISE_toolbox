function [y,newp,retcode]=sstate_model(obj,y,p,d,id) %#ok<INUSL>
% sstate_model -- shows the new way of writing a RISE steady state file
%
% Syntax
% -------
% ::
%
%   [y,newp,retcode]=sstate_model(obj,y,p,d,id)
%
% Inputs
% -------
%
% - **obj** [rise|dsge]: model object (not always needed)
%
% - **y** [vector]: endo_nbr x 1 vector of initial steady state
%
% - **p** [struct]: parameter structure
%
% - **d** [struct]: definitions
%
% - **id** [vector]: location of the variables to calculate
%
% Outputs
% --------
%
% - **y** []: endo_nbr x 1 vector of updated steady state
%
% - **newp** [struct]: structure containing updated parameters if any
%
% - **retcode** [0|number]: return 0 if there are no problems, else return
%   any number different from 0
%
% More About
% ------------
%
% - this is new approach has three main advantages relative to the previous
%   one:
%   - The file is valid whether we have many regimes or not
%   - The user does not need to know what regime is being computed
%   - It is in sync with the steady state model
%
% Examples
% ---------
%
% See also:

retcode=0;
if nargin==1
    
    y={'LG_C','LG_I','LG_H','Z_C','Z_I','ZL_C','ZL_I','ZG_C','ZG_I','AL',...
        'AG','A','G_C','G_I','G_H','H_C','K_C','K_I','I_C','I_I','I','H_I',...
        'H','C','LAMBDA_C','LAMBDA_I','XI_C','XI_I'};
    % flags on the calculation
    %--------------------------
    newp=struct('unique',true,'imposed',true,'initial_guess',false);
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