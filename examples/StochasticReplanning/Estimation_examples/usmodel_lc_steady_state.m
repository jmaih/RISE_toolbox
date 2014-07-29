function [ss,obj,retcode,imposed]=usmodel_lc_steady_state(obj,flag)

imposed=true;
retcode=0;
ss=[];
if flag==0
    params={'dy','dc','dinve','dw','pinfobs','robs','labobs'};
else
    pp=get(obj,'parameters');
	
	% the parameter below is not assigned a value in the model file
	% we set it and push it back into the rise object.
	nonset=struct('pelin_tp_2_1',1-pp.pelin_tp_1_2);
	obj=set(obj,'parameters',nonset);
    %---------------------------------
    cpie=1+pp.constepinf/100;
    cgamma=1+pp.ctrend/100 ;
    cbeta=1/(1+pp.constebeta/100);
    cr=cpie/(cbeta*cgamma^(-pp.csigma));
    conster=(cr-1)*100;
    %---------------------------------
    dy=pp.ctrend;
    dc=pp.ctrend;
    dinve=pp.ctrend;
    dw=pp.ctrend;
    pinfobs=pp.constepinf;
    robs = conster;
    labobs = pp.constelab;
    
    ss=[dy,dc,dinve,dw,pinfobs,robs,labobs]';
end
