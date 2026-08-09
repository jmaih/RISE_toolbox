%% housekeeping
clear
close all
clc

%% inspect the new model

% edit('occbin.rz')

%% rise the model

m = rise('occbin');

%% collect parameter values

p=baselineParameters();

p.r_zlb=1;

p.bind_ocb_1=0;

p.bind_ocb_2=1;

p.ocb_tp_1_2=0;
p.ocb_tp_2_1=0;

%% push parameters + steady state file

m=set(m,'parameters',p,'sstate_file','ssfile_full');

%% Inform RISE that we are solving a piecewise-linear model

m = set(m,'solve_occbin',1,... % this is the ID of the reference regime
    'sstate_imposed',true,... % the steady state is the steady state of the reference regime
    'sstate_unique',true); % the steady state is unique

%% solve the model

% N.B. the piecewise-linear machinery is assembled by the default solver.
% Do not override it here with e.g. 'solver',{'mfi',...}

m = solve(m);

%% print solution

print_solution(m)

%% the occasionally binding constraint

% Constraints are declared here and handed to simulate and irf, not written
% in the model file. The first column is the inequality that has to hold,
% the second the equation that holds while the constraint binds.

C={'R>=1','R=r_zlb'};

%%
clc

mysims = simulate(m,'simul_periods',1000,...
    'simul_constraints',C,...
    'simulator','occbin');

fprintf('lowest interest rate over the simulation: %0.10f\n',min(double(mysims.R)));

%% pick variables

var_list={'C','M','MC','N','PI','R','RR','W','Y'};

plot_sims(m,var_list,[],mysims)

%% IRF analysis

bigshock=5;

M=[set(m,'irf_shock_sign',-bigshock),...
    set(m,'irf_shock_sign',bigshock)];

myirfs=irf(M,'simul_constraints',C,'simulator','occbin');

%% plot each shock in turn

% N.B. impulse responses are not dated, so no date range is passed to
% plot_sims below.

shock_names=get(m,'exo_list');

tex=get(m,'tex');

for ii=1:numel(shock_names)

    shock=shock_names{ii};

    plot_sims(m,var_list,[],myirfs.(shock))

    [~,h]=sup_label(['Negative and positive ',tex.(shock)],'t');

    set(h,'fontsize',12)

    pause

    close

end
