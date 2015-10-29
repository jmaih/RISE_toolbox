function [eqtn,dico,leadorlag]=straigthen_equation(eqtn,variable,leadorlag,dico)

if abs(leadorlag)>0
    endo_names={dico.endogenous.name};
    param_exo=[{dico.parameters.name},{dico.exogenous.name}];
    
    is_param_Or_exo=any(strcmp(variable,param_exo));
    
    if is_param_Or_exo
        % create an auxiliary variable first and replace the location
        new_variable=parser.create_auxiliary_name(variable,0,true);
        % replace the item with the auxiliary variable
        eqtn{1,end}=new_variable;
        if ~any(strcmp(new_variable,endo_names))
            new_var=parser.listing('name',new_variable);
            %new_var=parser.listing('name',new_variable,'current_name',variable);
            dico.endogenous(end+1)=new_var;
            % correspondingly, an equation has to be created
            %------------------------------------------------
            create_latent_equation(new_variable,variable);
        end
        [eqtn,dico,leadorlag]=parser.straigthen_equation(eqtn,new_variable,leadorlag,dico);
        return
    elseif abs(leadorlag)>1
        % realm of endogenous variables
        vloc= strcmp(variable,endo_names);
        is_log_var=dico.endogenous(vloc).is_log_var;
        leadorlag__=leadorlag;
        dico.endogenous=parser.update_variable_lead_lag(dico.endogenous,...
            variable,leadorlag,is_log_var);
        if leadorlag>0
            new_var=parser.create_auxiliary_name(variable,leadorlag__-1,true);
            leadorlag=1;
        else
            new_var=parser.create_auxiliary_name(variable,leadorlag__+1,true);
            leadorlag=-1;
        end
        eqtn{1,end}=new_var;
        dico.endogenous=parser.update_variable_lead_lag(dico.endogenous,...
            new_var,leadorlag,is_log_var,variable);
    else
        is_log_var=false;
        vloc= strcmp(variable,endo_names);
        if any(vloc)
            is_log_var=dico.endogenous(vloc).is_log_var;
        end
        dico.endogenous=parser.update_variable_lead_lag(dico.endogenous,...
            variable,is_log_var,leadorlag);
    end
    eqtn{2,end}=leadorlag;
end

    function create_latent_equation(newguy,oldguy)
        newthing=sprintf('%s = %s;',newguy,oldguy);
        if isfield(dico,'latent_equations')
            dico.latent_equations=[dico.latent_equations
                {nan,newthing,'auxiliary equations'}];
        else
            dico.latent_equations={nan,newthing,'auxiliary equations'};
        end
    end

end