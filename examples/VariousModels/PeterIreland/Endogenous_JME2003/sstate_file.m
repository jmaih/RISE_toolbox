function [ys,obj,retcode,imposed]=sstate_file(obj,flag)

% This return code says that there is no problem so far
%------------------------------------------------------
retcode=0;

% This flag says that we do not impose the values calculated below. RISE
% will check whether they do solve the steady state
%--------------------------------------------------------------------------
imposed=false;

if flag==0
    % list the variables for which we compute the steady state. Not all
    % variables need to be listed. The variables NOT in this list
    % automatically get a value of 0. RISE is responsible for calculating
    % the steady state of any variable it may have created. The user does
    % not have to worry about those.
    %----------------------------------------------------------------------
    ys={'X','A','E','Z','V','MU','PAI','R','Q','LAMBDA','XI','C','M','Y',...
        'K','I','H','W','D','N','LC','LI','LM','LPI','LR'};
else
    % compute the steady state
    %-------------------------
    
    % get the parameters
    %--------------------
    p=get(obj,'parameters');
    
    % create definitions (auxiliary parameters) that are useful for solving
    % the steady state
    %----------------------------------------------------------------------
    z_ss=p.z_ss_trans*10000;
% 	phi_p = 100*abs(p.phi_p_trans);
% 	phi_k = 100*abs(p.phi_k_trans);
	mu_ss = 1 + p.mu_ss_trans;

    X=1;
    A=1;
    E=p.e_ss;
    Z=z_ss;
    V=1;
    % mu determined by policy
    MU=mu_ss;
    PAI=MU/p.g;
    R=MU/p.beta;
    Q=p.g/p.beta-1+p.delta;
    % lambda comes here
    l1ss=1/(p.theta/(p.theta-1)-(p.g-1+p.delta)*p.alpha/Q);
    l2ss=1/(1+E*(R/(R-1))^(p.gam-1));
    l3ss=(1-p.alpha)*Z*((p.theta-1)/p.theta)^(1/(1-p.alpha))*(p.alpha/Q)^(p.alpha/(1-p.alpha));
    LAMBDA=(p.eta+(1-p.alpha)*l1ss*l2ss)/l3ss;
    XI = (p.theta-1)/p.theta*LAMBDA;
    C = 1/(1+E*(R/(R-1))^(p.gam-1))*(1/LAMBDA);
    M = E*(R/(R-1))^p.gam*C;
    Y = 1/(1-(p.g-1+p.delta)*(p.theta-1)/p.theta*p.alpha/Q)*C;
    K = (p.theta-1)/p.theta*p.alpha*Y/Q;
    I = (p.g-1+p.delta)*K;
    H=1/Z*(Y/K^p.alpha)^(1/(1-p.alpha));
    W=(1-p.alpha)*(p.theta-1)/p.theta*Y/H;
    D = Y-W*H-Q*K;
    N=Y/H;
    LC=log(C);
    LI=log(I);
    LM=log(M);
    LPI=log(PAI);
    LR=log(R);
    ys=[X A E Z V MU PAI R Q LAMBDA XI C M Y K I H W D N LC LI LM LPI LR].';
    
    % if there is any problem in computing the steady state, let RISE know
    % about it. This is a flexibility we do not have if we write the steady
    % state inside the model file.
    %----------------------------------------------------------------------
    if any(ys>1e+5)||any(isnan(ys))||any(isinf(ys))
        retcode=1;
    end
end

end