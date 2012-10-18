function [ss,retcode]=steady_state_4_Canonical_Const(param_obj,flag)%param_struct=struct(parlist,parvals)
% if flag==0, then return the list of the variables, else return the values
% of their steady states

retcode=0;
switch flag
    case 0
        ss={'DPQ_P_NW','Y','ZGDP','ZI','ZPAI','ZY','D_GDP_NW','I','PAI','R','RN3M_NW'};
    case 1
        number_of_regimes=numel(param_obj(1).startval);
        number_of_variables=11;
        param_names={param_obj.name};
        
        par_vals=vertcat(param_obj.startval);
        ss=nan(number_of_variables,number_of_regimes);
        for ii=1:number_of_regimes
            [ss(:,ii)]=regime_specifi_steady_state(param_names,par_vals(:,ii));
        end
    otherwise
        error([mfilename,':: unknown flag'])
end


function ss=regime_specifi_steady_state(param_names,params)  %#ok<INUSD>
for index=1:numel(param_names)
	eval([param_names{index},'=params(index);'])
end

DPQ_P_NW=paiss;
Y=0;
ZGDP=gyss;
ZI=0;
ZPAI=0;
ZY=0;
D_GDP_NW=gyss;
I=0;
PAI=0;
R=0;
RN3M_NW=iss;

ss=[DPQ_P_NW , Y , ZGDP , ZI , ZPAI , ZY , D_GDP_NW , I , PAI , R , RN3M_NW ]' ;
