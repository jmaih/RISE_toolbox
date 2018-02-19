%% housekeeping
close all
clear
clc
%% load models

load tut02_estimation % tut_estim_combos

model_names=fieldnames(models);

%% set shocks

shock_names={'fmp','productivity','cost_push','mp','FX'};

%% set up restrictions: replicating the choleski restrictions

ident_restr1={
    % normalization restrictions
    'rtwi{0}@fmp','+'
    'LGDP{0}@productivity','+'
    'Aldp{0}@cost_push','+'
    'r{0}@mp','+'
    'LRER{0}@FX','+'
    % first set
    'rtwi{0}@productivity',0
    'rtwi{0}@cost_push',0
    'rtwi{0}@mp',0
    'rtwi{0}@FX',0
    % second set
    'LGDP{0}@cost_push',0
    'LGDP{0}@mp',0
    'LGDP{0}@FX',0
    % third set
    'Aldp{0}@mp',0
    'Aldp{0}@FX',0
    % fourth set
    'r{0}@FX',0 % This last restriction, we will vary
    };

%% Choleski (top) + Long run restriksjon på delta valutakurs

ident_restr2=ident_restr1;

ident_restr2(end,:)={'LRER{inf}@mp',0};

%% Choleski (top) + sign på valutakurs

ident_restr3=ident_restr1;

ident_restr3(end,:)={'LRER{0}@mp','-'};

%% store restrictions for later use
								
agnostic=true;

max_trials=6000;

Rfunc1=struct(); ident1=struct();

Rfunc2=struct(); ident2=struct();

Rfunc3=struct(); ident3=struct();

for imod=1:numel(model_names)
    
    m=model_names{imod};
    
    [Rfunc1.(m),ident1.(m)]=identification(models.(m),ident_restr1,shock_names,...
        agnostic,max_trials);
    
    [Rfunc2.(m),ident2.(m)]=identification(models.(m),ident_restr2,shock_names,...
        agnostic,max_trials);
    
    [Rfunc3.(m),ident3.(m)]=identification(models.(m),ident_restr3,shock_names,...
        agnostic,max_trials);

end

%% save for later use

save('tut04_identification','Rfunc1','Rfunc2','Rfunc3','ident1','ident2','ident3','shock_names')
