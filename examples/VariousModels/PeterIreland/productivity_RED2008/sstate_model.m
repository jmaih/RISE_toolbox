function [ys,obj,retcode,imposed]=sstate_model(obj,flag)

retcode=0;

imposed=false;

if flag==0
    
    ys={'LG_C','LG_I','LG_H','Z_C','Z_I','ZL_C','ZL_I','ZG_C','ZG_I','AL',...
        'AG','A','G_C','G_I','G_H','H_C','K_C','K_I','I_C','I_I','I','H_I',...
        'H','C','LAMBDA_C','LAMBDA_I','XI_C','XI_I'};
else
    p=get(obj,'parameters');
    
    defs=get(obj,'definitions');

    newp=struct();
    
    newp.theta_I=p.theta_C;
    
	AL=1;

	AG=defs.ag;
   
	A=AL*AG;

    G_C = AG*(defs.zg_I^p.theta_C)*(defs.zg_C^(1-p.theta_C));
    
    G_I = AG*defs.zg_I;
    
    G_H = AG;
    
    H_C = (1-p.theta_C)*((G_C-p.beta*p.gam)/(G_C-p.gam))*AG;
    
    kcssa = (1-newp.theta_I)*AG*defs.zg_I;
    
    kcssb = (G_C-p.beta*p.gam)/(G_C-p.gam);
    
    kcssc = p.theta_C*(AG*defs.zg_I-p.beta*(1-p.delta_I));
    
    kcssd = newp.theta_I*(AG*defs.zg_I-p.beta*(1-p.delta_C));
    
    kcsse = ((p.beta*newp.theta_I)/(AG*defs.zg_I-p.beta*(1-p.delta_I)))^(1/(1-newp.theta_I));
    
    K_C = kcssa*kcssb*(kcssc/kcssd)*kcsse;
    
    kissa = p.beta*newp.theta_I*(AG*defs.zg_I-1+p.delta_C);
    
    kissb = AG*defs.zg_I - p.beta*(1-p.delta_I);
    
    kissc = p.beta*newp.theta_I*(AG*defs.zg_I-1+p.delta_I);
    
    K_I = (kissa/(kissb-kissc))*K_C;
    
    I_C = (AG*defs.zg_I-1+p.delta_C)*K_C;
    
    I_I = (AG*defs.zg_I-1+p.delta_I)*K_I;
    
    I = I_C + I_I;
    
    H_I = (1/defs.zg_I)*((I/(K_I^newp.theta_I))^(1/(1-newp.theta_I)));
    
    H = H_C + H_I;
    
    C = (K_C^p.theta_C)*((defs.zg_C*H_C)^(1-p.theta_C));
    
    LAMBDA_C = ((G_C-p.beta*p.gam)/(G_C-p.gam))/C;
    
    LAMBDA_I = H_I/((1-newp.theta_I)*AG*I);
    
    XI_C = LAMBDA_I;
    
    XI_I = XI_C;
    
	ZL_C=1;

	ZL_I=1;

	ZG_C=defs.zg_C;
    
	ZG_I=defs.zg_I;
    
    Z_C=ZL_C*ZG_C;
	
    Z_I=ZL_I*ZG_I;
    
    LG_C=log(G_C);
    
    LG_I=log(G_I);
    
    LG_H=log(G_H);
	
    newp.eta_C=AG;
    
    newp.eta_I=AG; 
    
    newp.kappa_C=I_C/K_C; 
    
    newp.kappa_I=I_I/K_I;
    
    obj=set(obj,'parameters',newp);
    
    ys=[LG_C,LG_I,LG_H,Z_C,Z_I,ZL_C,ZL_I,ZG_C,ZG_I,AL,AG,A,G_C,G_I,G_H,H_C,...
        K_C,K_I,I_C,I_I,I,H_I,H,C,LAMBDA_C,LAMBDA_I,XI_C,XI_I];
    
    ys=ys(:);
    
    retcode=any(isnan(ys))||any(isinf(ys));
end

end