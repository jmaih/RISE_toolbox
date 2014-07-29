% computes the steady state for the observed variables in the smets-wouters
% model. You only need to provide the steady state for the variables whose
% steady state is different from zero.
function [ys,retcode]=usmodel_steadystate(param_obj,flag)

retcode=0;
switch flag
    case 0
        ys={'dy','dc','dinve','dw','pinfobs','robs','labobs'};
    case 1
        param_names={param_obj.name};
        params=vertcat(param_obj.startval); %#ok<NASGU>
        
        for index=1:numel(param_names)
            eval([param_names{index},'=params(index);'])
        end
        
        % In the SW model, one of the steady state is endogenous...
        cpie=1+constepinf/100;
        cgamma=1+ctrend/100 ;
        cbeta=1/(1+constebeta/100);
        %         clandap=cfc;
        %         cbetabar=cbeta*cgamma^(-csigma);
        cr=cpie/(cbeta*cgamma^(-csigma));
        %         crk=(cbeta^(-1))*(cgamma^csigma) - (1-ctou);
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
        
        
        dy=ctrend;
        dc=ctrend;
        dinve=ctrend;
        dw=ctrend;
        pinfobs = constepinf;
        robs =conster;
        labobs =constelab;
        
        ys =[dy,dc,dinve,dw,pinfobs,robs,labobs]';
    otherwise
        error([mfilename,':: Unknown flag'])
end

