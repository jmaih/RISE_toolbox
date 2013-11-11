function [ys,obj,retcode,imposed]=usmodel_steadystate(obj,flag)

% computes the steady state for the observed variables in the smets-wouters
% model. You only need to provide the steady state for the variables whose
% steady state is different from zero.

imposed=true;
retcode=0;
switch flag
    case 0
        ys={'dy','dc','dinve','dw','pinfobs','robs','labobs'};
    case 1
        pp=get(obj,'parameters');
        
        % In the SW model, one of the steady state is endogenous...
        cpie=1+pp.constepinf/100;
        cgamma=1+pp.ctrend/100 ;
        cbeta=1/(1+pp.constebeta/100);
        %         clandap=cfc;
        %         cbetabar=cbeta*cgamma^(-pp.csigma);
        cr=cpie/(cbeta*cgamma^(-pp.csigma));
        %         crk=(cbeta^(-1))*(cgamma^pp.csigma) - (1-ctou);
        %         cw = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*crk^calfa))^(1/(1-calfa));
        %         cikbar=(1-(1-ctou)/cgamma);
        %         cik=(1-(1-ctou)/cgamma)*cgamma;
        %         clk=((1-calfa)/calfa)*(crk/cw);
        %         cky=cfc*(clk)^(calfa-1);
        %         ciy=cik*cky;
        %         ccy=1-cg-cik*cky;
        %         crkky=crk*cky;
        %         cwhlc=(1/clandaw)*(1-calfa)/calfa*crk*cky/ccy;
        %         cwly=1-crk*cky;
        conster=(cr-1)*100;
        
        
        dy=pp.ctrend;
        dc=pp.ctrend;
        dinve=pp.ctrend;
        dw=pp.ctrend;
        pinfobs = pp.constepinf;
        robs =conster;
        labobs =pp.constelab;
        
        ys =[dy,dc,dinve,dw,pinfobs,robs,labobs]';
        if any(isinf(ys))||any(isnan(ys))
            % if there is a problem, return a flag different from zero. It
            % does not matter what number you write as long as it is
            % different from zero.
            retcode=1;
        end
    otherwise
        error([mfilename,':: Unknown flag'])
end