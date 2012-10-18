function [ss,var_names,retcode]=steady_state_4_Canonical_Const(param_obj)%param_struct=struct(parlist,parvals)

persistent param_names number_of_regimes number_of_variables

retcode=0;
if isempty(number_of_regimes)
number_of_regimes=numel(param_obj(1).startval);
	number_of_variables=11;
	param_names={param_obj.name};
end

par_vals=vertcat(param_obj.startval);
ss=nan(number_of_variables,number_of_regimes);
var_names=[];
for ii=1:number_of_regimes
	[ss(:,ii),var_names]=regime_specifi_steady_state(param_names,par_vals(:,ii),var_names);
end


function [ss,var_names]=regime_specifi_steady_state(param_names,params,var_names) %#ok<INUSL>
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
if isempty(var_names)
    var_names={'DPQ_P_NW','Y','ZGDP','ZI','ZPAI','ZY','D_GDP_NW','I','PAI','R','RN3M_NW'};
end