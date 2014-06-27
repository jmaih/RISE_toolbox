%% housekeeping
clear all
close all
clc
%% rise the model file

gdexo=rise('globdynzlb','rise_flags',struct('exogenous_switching',true));
gdexo.legend='exo_switch';

% gdendo=rise('globdynzlb','rise_flags',struct('exogenous_switching',false));
% gdendo.legend='endo_switch';
% % 
% gd=[gdexo,gdendo];
%% solve the model

gdexo=solve(gdexo);

%% print the solution

print_solution(gdexo)

%% compute irfs

myirfs=irf(gdexo,'irf_risk',0,'irf_periods',20);
mygirfs=irf(gdexo,'irf_type','girf','irf_periods',10,'irf_risk',0);
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
        if ishock==1 && ivar==1
            nobs=myirfs.(shock).(vname).NumberOfObservations;
            ssbench=rise_time_series.ones(0,nobs,1);
        end
        hold on
        plot(full(gdexo.solution.ss{1}(vlocs(ivar)))*ssbench,...
            'linewidth',2,'color',[1 0 0])
        title([vtexNames{ivar},' (on ',shock_texnames{ishock},')'])
    end
end

%% the system is mean-square stable
for ishock=1:numel(shock_names)
    shock=shock_names{ishock};
    figure('name',['Generalized impulse responses to a ',shock_texnames{ishock}])
    for ivar=1:numel(varlist)
        vname=varlist{ivar};
        ssvar=full(gdexo.solution.ss{1}(vlocs(ivar)));  
        subplot(3,3,ivar)
        plot(ssvar+mygirfs.(shock).(vname),'linewidth',2),
        if ishock==1 && ivar==1
            nobs=mygirfs.(shock).(vname).NumberOfObservations;
            ssbench=rise_time_series.ones(0,nobs,1);
        end
        hold on
        plot(ssvar*ssbench,...
            'linewidth',2,'color',[1 0 0])
        title([vtexNames{ivar},' (on ',shock_texnames{ishock},')'])
    end
end
%% we can also solve the model using a 2nd-order perturbation

ldzlb2=solve(gdexo,'solve_order',2);

%% and print the solution

print_solution(ldzlb2)

%% select the variables to print
print_solution(ldzlb2,{'C','I','K','N','PAI','PSI','R','RK'})

%% print the solution in a more compact form
print_solution(ldzlb2,{'C','I','K','N','PAI','PSI','R','RK'},true)

