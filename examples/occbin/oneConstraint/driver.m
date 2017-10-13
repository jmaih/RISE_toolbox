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

m=set(m,'parameters',p,'steady_state_file','ssfile_full');

%% Inform RISE that we are solving a piecewise-linear model

m = set(m,'solve_occbin',1,... % this is the ID of the reference regime
    'steady_state_imposed',true,... % the steady state is the steady state of the reference regime
    'steady_state_unique',true); % the steady state is unique

%% solve the model

m = solve(m);

%% print solution

print_solution(m)

%%
clc

mysims = simulate(m,'simul_periods',1000,...
    'simul_honor_constraints',true);

%% pick variables

var_list={'C','M','MC','N','PI','R','RR','W','Y'};

plot_sims(m,var_list,[],mysims)

%% IRF analysis
bigshock=5;

M=[set(m,'irf_shock_sign',-bigshock),...
    set(m,'irf_shock_sign',bigshock)];

%% do the irfs without honoring constraints

myirfs_nc=irf(M);

%% plot each shock in turn
shock_names=get(m,'exo_list');

tex=get(m,'tex');

for ii=1:numel(shock_names)
    
    shock=shock_names{ii};
    
    plot_sims(m,var_list,'0:15',myirfs_nc.(shock))
    
    [~,h]=sup_label(['Negative and positive ',tex.(shock)],'t');
    
    set(h,'fontsize',12)
    
    pause
    
    close
    
end

%% do the irfs honoring constraints

myirfs_c=irf(M,'simul_honor_constraints',true);

%% plot each shock in turn
shock_names=get(m,'exo_list');

tex=get(m,'tex');

for ii=1:numel(shock_names)
    
    shock=shock_names{ii};
    
    plot_sims(m,var_list,'0:15',myirfs_c.(shock))
    
    [~,h]=sup_label(['Negative and positive ',tex.(shock)],'t');
    
    set(h,'fontsize',12)
    
    pause
    
    close
    
end

