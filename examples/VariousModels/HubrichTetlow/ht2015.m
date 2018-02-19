%% housekeeping
clear
clc
close all

%% create data

tmp=importdata('ht2014data.txt');

vnames={'real_pce','cpixfe','ffr','m2','fsi'};

rawdb=ts('1988M12',tmp.data(:,2:end),vnames);

rawdb=pages2struct(rawdb);

%% plot the data

figure('name','raw data');
for ii=1:numel(vnames)
    v=vnames{ii};
    subplot(3,2,ii)
    plot(rawdb.(v))
    title(v)
end

%% transformed data
% 
% we use levels of the federal funds rate and the stress index and growth
% rates of real personal consumption expenditures (PCE), money and prices.
db=struct();
db.R=rawdb.ffr;
db.S=rawdb.fsi;
db.C=100*log(rawdb.real_pce/rawdb.real_pce{-1});
db.M=100*log(rawdb.m2/rawdb.m2{-1});
db.P=100*log(rawdb.cpixfe/rawdb.cpixfe{-1});

descr=struct('C','Consumption growth',...
    'P','Inflation',...
    'R','Feds Funds rate',...
    'M','Money growth',...
    'S','Financial stress index');

varlist=fieldnames(descr);

figure('name','Transformed data');
for ii=1:numel(varlist)
    v=varlist{ii};
    subplot(3,2,ii)
    plot(db.(v))
    title(descr.(v))
    axis tight
end
%% Model 1: SVAR models with cholesky identification
clc

models={'1v1c','2v1c','3v1c','1v2c','2v2c','3v2c',...
    '3vS2c','3vSC2c','3vSCP2c','3vSRM2c','3vRM2c','3vRMC2c'};

nlags=2;

exog={};

constant=true;

panel=[];

priors=struct();
priors.var=svar.prior_template();
priors.var.type='sz';

date_range={db.C.start,db.C.finish};

sv0=svar.empty(1,0);

tic

for ii=1:numel(models)
    
    fprintf(1,' -------------------- %s -------------------- \n',models{ii});
    
    [markov_chains,switch_prior,restrictions]=build_model(models{ii},varlist);
    
    sv0(1,ii)=svar(varlist,exog,nlags,constant,panel,markov_chains);
    
    priors.nonvar=switch_prior;
    
    sv0(1,ii)=estimate(sv0(1,ii),db,date_range,priors,restrictions);
    
end

toc

% 2601.732325 seconds

%% estimates

for ii=1:numel(models)
    clc, close all
    
    fprintf(1,' -------------------- %s -------------------- \n',models{ii});
    
    pmode=posterior_mode(sv0(ii))
    
    pause
    % Printing estimates
    
    print_structural_form(sv0(ii))
    
    pause
    % Printing solution
    
    print_solution(sv0(ii))
    pause

    % historical probabilities
    plot_probabilities(sv0(ii))
    pause
    
    % plots probabilities against data
    close all
    plot_data_against_probabilities(sv0(ii),'regime')
    pause
    
end

%% A model with endogenous switching

clc

[mctvp,switch_prior,restrictions]=build_model('2v1c',varlist);

mctvp.endogenous_probabilities={
    'syncvol_tp_1_2=1/(1+exp(a12-b12*S))'
    'syncvol_tp_2_1=a21'
    };

mctvp.probability_parameters={'a12','b12','a21'};

mctvp.controlled_parameters={'s(3)'};%'s(1)','s(2)','s(4)',

switch_prior=rmfield(switch_prior,{'syncvol_tp_1_2','syncvol_tp_2_1'});

switch_prior.a12={0.1,0,5,'normal'};

switch_prior.b12={0.1,2,5,'gamma'};

switch_prior.a21={0.2,0.2,0.2,'beta'};

svtvp=svar(varlist,exog,nlags,constant,panel,mctvp);
    
priors.nonvar=switch_prior;
    
svtvp=estimate(svtvp,db,date_range,priors,restrictions);

