%% Housekeeping
clear all
close all
clc
%% load RISE
rise_startup()


%% read the models and their calibrations
sw=rise('usmodel');

%% solve the model
sw=solve(sw,'solver','newton_kronecker_iteration');

%% print results
sw.print_solution()

%% print solution for a subset of variables only
sw.print_solution({'a','b','c','cf','dc','dinve','dw','dy'})

%% compute regime-specific impulse responses
myirfs0=irf(sw,'irf_periods',20);

%% compute generalized impulse responses
myirfs1=irf(sw,'irf_periods',20,'irf_type','girf');

%% plot the impulse responses
close all
var_list={'dc','dinve','dw','dy','lab','pinf','r'};
figure('name','Impulse responses to a wage markup shock');
for ii=1:numel(var_list)
    subplot(3,3,ii)
	v=var_list{ii};
	aggregate=[myirfs0.ew.(v),myirfs1.ew.(v)];
    plot(aggregate,'linewidth',2);
    title(var_list{ii})
    if ii==1
        legend(aggregate.varnames)
    end
    axis tight
end

