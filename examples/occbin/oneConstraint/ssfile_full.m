function [y,newp,retcode]=ssfile_full(obj,y,p,d,id)

retcode=0;

if nargin==1
    
    y={'PI','MUZ','R_TAYLOR','A','Y','C','LAMBDA','PIND','M','MC','N',...
        'W','DLOG_Y','RR','R'};
    
    newp='kappa';
    
else
    
    newp=struct();
    
    PI = p.pi_ss;
    
    MUZ = exp(p.gz);
    
    R_TAYLOR = p.pi_ss*exp(p.gz)/p.beta;
    
    R = R_TAYLOR;
    
    A = 1;
    
    Y = 1;
    
    C = 1;
    
    LAMBDA = (C - p.h*C/MUZ)^(-p.sigma);
    
    PIND = PI;
    
    M = 1/R;
    
    MC = (p.epsilon - 1)/p.epsilon;
    
    N = 1;
    
    W = MC;
    
    kappa = LAMBDA*W/(N^p.eta);
    
    DLOG_Y = log(MUZ);
    
    RR = R/PI;
    
    ytemp=[PI,MUZ,R_TAYLOR,A,Y,C,LAMBDA,PIND,M,MC,N,W,DLOG_Y,RR,R];
    
    y(id)=ytemp;
    
    newp.kappa=kappa;
    
    if any(isnan(y)|isinf(y))
        
        retcode=1;
        
    end
    
end


end
