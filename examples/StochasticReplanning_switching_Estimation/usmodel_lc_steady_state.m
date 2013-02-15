function [params,ss,retcode,imposed]=usmodel_lc_steady_state(params,flag)

imposed=true;
retcode=0;
ss=[];
if flag==0
    params={'dy','dc','dinve','dw','pinfobs','robs','labobs'};
else
    name_loc=strcmp('name',params(1,:));
    val_loc=strcmp('startval',params(1,:));
    par_names=params{2,name_loc};
    par_mat=params{2,val_loc};
    parlist={'ctrend','constepinf','conster','constebeta','constelab',...
        'csigma','pelin_tp_1_2','pelin_tp_2_1'};
    pp=struct();
    for ilist=1:numel(parlist)
        par=parlist{ilist};
        loc=strcmp(par,par_names);
        if strcmp(par,'pelin_tp_2_1')
            % push back the parameter inside the params
            params{2,val_loc}(loc,:)=1-pp.pelin_tp_1_2;
            continue
        end
        pp.(par)=par_mat(loc,1);
    end
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
