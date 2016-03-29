%% housekeeping
close all
clc
%% rise the model file

gdexo=rise('globdynzlb','rise_flags',struct('exogenous_switching',true));
gdexo.legend='endo_switch';

%% solve the model

gdexo=solve(gdexo,'steady_state_imposed',true);

%% print the solution

print_solution(gdexo)

%% compute irfs

myirfs=irf(gdexo,'irf_periods',20);
mygirfs=irf(gdexo,'irf_type','girf','irf_periods',10,'irf_regime_specific',false);
%% plot irfs
endo_names=gdexo.endogenous.name;
endo_texnames=gdexo.endogenous.tex_name;
shock_names=gdexo.exogenous.name;
shock_texnames=gdexo.exogenous.tex_name;
varlist={'C','PAI','R','N','Y','R','W','PSI','PSI'};
vlocs=locate_variables(varlist,gdexo.endogenous.name);
vtexNames=gdexo.endogenous.tex_name(vlocs);

for ishock=1:numel(shock_names)
    shock=shock_names{ishock};
    figure('name',['impulse responses to a ',shock_texnames{ishock}])
    for ivar=1:numel(varlist)
        vname=varlist{ivar};
        subplot(3,3,ivar)
        plot(myirfs.(shock).(vname),'linewidth',2)
        if ivar==1
            legend('normal','zlb')
        end
        title([vtexNames{ivar},' (on ',shock_texnames{ishock},')'])
    end
end

%% the system is mean-square stable
for ishock=1:numel(shock_names)
    shock=shock_names{ishock};
    figure('name',['Generalized impulse responses to a ',shock_texnames{ishock}])
    for ivar=1:numel(varlist)
        vname=varlist{ivar};
        subplot(3,3,ivar)
        plot(mygirfs.(shock).(vname),'linewidth',2),
        if ishock==1 && ivar==1
            nobs=mygirfs.(shock).(vname).NumberOfObservations;
        end
    end
end
%% we can also solve the model using a 2nd-order perturbation

ldzlb2=solve(gdexo,'solve_order',2,'solve_derivatives_type','automatic');

%% and print the solution

print_solution(ldzlb2)

%% select the variables to print
print_solution(ldzlb2,{'C','I','K','N','PAI','PSI','R','RK'})
