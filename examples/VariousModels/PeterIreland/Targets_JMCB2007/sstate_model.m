function [ys,obj,retcode,imposed]=sstate_model(obj,flag)

imposed=false;
retcode=0;
if flag==0
    ys={'A','THETA','Z','PAISTAR','V','PAI','Y','C','GY','GPAI','GR',...
        'LAMBDA','Q','X','R','RRPAI','E'};
else
    p=get(obj,'parameters');
    A=1;
    THETA=p.thetass;
    Z=p.zss;
    PAISTAR=1;
    V=1;
    PAI=1;
    Y=(p.thetass-1)/p.thetass*(p.zss-p.beta*p.gam)/(p.zss-p.gam);
    C=Y;
    GY=Z;
    GPAI=1;
    GR=1;
    LAMBDA=p.thetass/(p.thetass-1);
    Q=(p.zss-p.beta*p.gam)/(p.zss-p.gam);
    X=(p.thetass-1)/p.thetass;
    R=p.zss/p.beta;
    RRPAI=p.zss/p.beta;
    E=p.ess;
    
    ys=[A,THETA,Z,PAISTAR,V,PAI,Y,C,GY,GPAI,GR,LAMBDA,Q,X,R,RRPAI,E];
    ys=ys(:);
    if any(isnan(ys))||any(isinf(ys))
        retcode=1;
    end
end

end