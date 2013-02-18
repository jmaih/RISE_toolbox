%% housekeeping
clear all
close all
clc

%% read the model file

lwz=rise('lwz09_2','irf_periods',20);

%% the naive economy: assume the Hawkish regime lasts forever
% q_tp_2_1=0 for all regimes
naive_calibration=struct();
naive_calibration.q_tp_2_1=0;

lwz_naive=lwz.set_parameters(naive_calibration);
%% solve, evaluate, print solution, etc, not required

%% compute impulse responses
myirfs=irf(lwz);

%% compute irfs for the naive economy
my_naive_irfs=irf(lwz_naive);
%% plot impulse responses: naive vs normal

shock_names={lwz.varexo.name};
mylist={'PAI','Y','R'};
for ishock=1:numel(shock_names)
    shock=shock_names{ishock};
    figure('name',['IRFs to a ',lwz.varexo(ishock).name,' shock (',shock,')'])
    for ivar=1:numel(mylist)
        vname=mylist{ivar};
        subplot(3,1,ivar)
        expect=myirfs.(shock).(vname);
        naive=my_naive_irfs.(shock).(vname);
        tmp=plot([expect,naive],'linewidth',2);
        set(tmp(1:2),'linestyle','--')
        title(vname)
        if ivar==1
            legend({'dovish','hawkish','naive-dovish','naive-hawkish'})
        end
    end
end

%% regime-dependent structural parameters ( eta and iota)

regDep_calibration=struct();
regDep_calibration.eta_q_2=0.75; % eta controled by q assumes 0.75 in state 2
regDep_calibration.iota_q_2=0;

lwz=lwz.set_parameters(regDep_calibration);
lwz_naive=lwz_naive.set_parameters(regDep_calibration);

%% impulse responses
myirfs_2=irf(lwz);

my_naive_irfs_2=irf(lwz_naive);

%% plot IRFs
close all
for ishock=1:numel(shock_names)
    shock=shock_names{ishock};
    figure('name',['IRFs (Regime dep. structural params) to a ',lwz.varexo(ishock).name,' shock (',shock,')'])
    for ivar=1:numel(mylist)
        vname=mylist{ivar};
        subplot(3,1,ivar)
        expect=myirfs.(shock).(vname);
        naive=my_naive_irfs.(shock).(vname);
        tmp=plot([expect,naive],'linewidth',2);
        set(tmp(1:2),'linestyle','--')
        title(vname)
        if ivar==1
            legend({'dovish','hawkish','naive-dovish','naive-hawkish'})
        end
    end
end
