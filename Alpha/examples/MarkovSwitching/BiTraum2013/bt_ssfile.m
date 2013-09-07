function [ss,obj,retcode,imposed]=bt_ssfile(obj,flag)

retcode=0;
% Instruct RISE not to bother checking that this is the true steady state
%------------------------------------------------------------------------
% alternatively, RISE will use the values computed herein as a starting
% value
imposed=true;
if flag==0
    % when flag is 0, return the list of the variables for which the steady
    % state will be computed
    ss={'N','A','Y','TAU','C','G','Q','B','Z','S','BD'};
else
    % when flag is not 0, compute the steady state
    
    % we get the parameters in a structure
    %-------------------------------------
    pp=get(obj,'par_vals');
    
    delta_ss=0;
    % we use the parameters to compute steady state values
    %-----------------------------------------------------
    N=pp.n;
    A=pp.a;
    Y=A.*N;
    TAU=pp.tau;
    G=pp.g_y.*Y;
    C=Y-G;
    B=pp.b_y.*Y;
    Z=(TAU-B./Y.*(1-pp.beta).*(1-delta_ss)-G./Y).*Y;
    Q=pp.beta.*(1-delta_ss);
    S=B./Y;
    BD=(1-delta_ss).*B;
    ss=[
        N
        A
        Y
        TAU
        C
        G
        Q
        B
        Z
        S
        BD
        ];
    
    % Let RISE know of any problem resulting in an unsuccessful calculation
    % of the steady state. So it does not proceed to solving trying to
    % solve the model
    %----------------------------------------------------------------------
    if any(~isfinite(ss(:)))
        retcode=1;
    end
    
    if ~retcode
        % some parameters determined as a function of the steady state
        %-------------------------------------------------------------
        shat = pp.stilde - 0.6*Y;
        
        eta2 = 1./(pp.stilde-shat).*log(pp.ptilde.*(1-pp.phat)./(pp.phat.*(1-pp.ptilde))); %  (15) p.11
        
        eta1 = log(pp.ptilde./(1-pp.ptilde)) - eta2.*pp.stilde; %  (15) p.11
        
        % they need to be pushed back into the model. Else solving may fail
        %------------------------------------------------------------------
        % we know that all those parameters are controlled by the constant
        % markov chain.
        pout=struct('shat',shat(1),'eta2',eta2(1),'eta1',eta1(1));
        obj=set(obj,'parameters',pout);
    end
end