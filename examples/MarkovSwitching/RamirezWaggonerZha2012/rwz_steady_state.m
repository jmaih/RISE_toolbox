function [ss,retcode,imposed]=rwz_steady_state(param_obj,flag)%param_struct=struct(parlist,parvals)
% if flag==0, then return the list of the variables, else return the values
% of their steady states

imposed=true;

retcode=0;
switch flag
    case 0
        ss={'PAI','Y','R'};
    case 1
        number_of_regimes=numel(param_obj(1).startval);
        number_of_variables=3;
        param_names={param_obj.name};
        
        par_vals=vertcat(param_obj.startval);
        ss=nan(number_of_variables,number_of_regimes);
        for ii=1:number_of_regimes
            [ss(:,ii)]=regime_specific_steady_state(param_names,par_vals(:,ii));
        end
    otherwise
        error([mfilename,':: unknown flag'])
end


function ss=regime_specific_steady_state(param_names,params)  %#ok<INUSD>
for index=1:numel(param_names)
	eval([param_names{index},'=params(index);'])
end

PAI=1;
Y=(eta-1)/eta;
R=exp(mu_bar)/betta*PAI;

ss=[PAI,Y,R]' ;
