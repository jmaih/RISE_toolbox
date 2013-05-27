function [ss,param_obj,retcode,imposed]=steady_state_4_Canonical_Const(param_obj,flag)%param_struct=struct(parlist,parvals)
% if flag==0, then return the list of the variables, else return the values
% of their steady states

% below, we lay out one way of computing the steady state with different
% parameterizations in e.g. the case where there are multiple regimes with
% specific steady states.

retcode=0;
imposed=false;
switch flag
    case 0
        ss={'DPQ_P_NW','Y','ZGDP','ZI','ZPAI','ZY','D_GDP_NW','I','PAI','R','RN3M_NW'};
    case 1
        pp=struct();
        name_loc=strcmp('name',param_obj(1,:));
        val_loc=strcmp('startval',param_obj(1,:));
        par_names=param_obj{2,name_loc};
        par_mat=param_obj{2,val_loc};
        number_of_regimes=size(par_mat,2);
        number_of_variables=11;
        ss=nan(number_of_variables,number_of_regimes);
        for ireg=1:number_of_regimes
            for ipar=1:numel(par_names)
                pp(ireg).(par_names{ipar})=par_mat(ipar,1);
            end
            [ss(:,ireg)]=regime_specifi_steady_state(pp(ireg));
        end
    otherwise
        error([mfilename,':: unknown flag'])
end


function ss=regime_specifi_steady_state(pp)  

DPQ_P_NW=pp.paiss;
Y=0;
ZGDP=pp.gyss;
ZI=0;
ZPAI=0;
ZY=0;
D_GDP_NW=pp.gyss;
I=0;
PAI=0;
R=0;
RN3M_NW=pp.iss;

ss=[DPQ_P_NW , Y , ZGDP , ZI , ZPAI , ZY , D_GDP_NW , I , PAI , R , RN3M_NW ]' ;
