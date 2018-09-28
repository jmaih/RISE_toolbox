%% housekeeping
clear
close all
clc
%% bring in the data
[num,txt] = xlsread('Smets_Wouters_data.xlsx');
vnames=strrep(txt,' ','');
date_start='1947q3';
data=ts(date_start,num(:,2:end),vnames);
data=pages2struct(data);
figure('name','observed data')
for ii=1:numel(vnames)
    subplot(3,3,ii)
    v=vnames{ii};
    plot(data.(v),'linewidth',2)
    title(v)
end
xrotate(45)

%% Set the RFVAR
endo_aliases=struct('robs','interest rates',...
    'dy','GDP growth',...
    'labobs','hours worked',...
    'pinfobs','inflation',...
    'dw','wages growth',...
    'dc','Consumption growth',...
    'dinve','Investment growth');

endog=fieldnames(endo_aliases);

nlags=3;

const=true;

exog={};

rv=rfvar(endog,exog,nlags,const); % formerly redfvar

%% Estimate the reduced-form VAR

rv1=estimate(rv,data);%,{db.LGDP.start,db.LGDP.finish}

%% identification setup

shock_aliases=struct('mp','monetary policy',...
    'ad','aggreg. demand',...
    'as','aggreg. supply',...
    'ls','labor supply');

structural_shocks=fieldnames(shock_aliases);

ident_restr={
    %     'dy{inf}@mp',0
    %     'dy{inf}@ad',0
    %     'dy{inf}@as',0
    %     'dy{inf}@ls',0
    'robs{0}@mp','+'
    'dy@mp','-'
    'pinfobs@mp','-'
    'robs{0}@ad','+'
    'pinfobs@ad','+'
    'dy@ad','+'
    %     'robs@as','-'
    'dy@as','+'
    'pinfobs@as','-'
    'robs@ls','-'
    'dy@ls','+'
    'labobs@ls','+'
    %     'pinfobs@ls','-'
    'dw@ls','-'
    };

agnostic=true;

max_trials=6000;

[Rfunc,ident]=identification(rv1,ident_restr,structural_shocks,...
        agnostic,max_trials);

%% Impulse responses

myirfs=irf(rv1,structural_shocks,40,[],Rfunc);

%% Plot the IRFs for one particular rotation

for ishock=1:numel(structural_shocks)
    
    shock=structural_shocks{ishock};
    
    figure('name',['Orthogonalized responses to a ',...
        shock_aliases.(shock),' shock']);
    
    for ivar=1:numel(endog)
        
        vname=endog{ivar};
        
        subplot(3,3,ivar)
        
        plot(myirfs.(shock).(vname),'linewidth',2)
        
        title(endo_aliases.(vname))
        
    end
    
end

