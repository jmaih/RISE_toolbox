%% -------------------- only variances switch model -------------------- %%
%
% Only variances in ALL three equations switch

%% housekeeping
clear
close all
clc()

%% Create dataset
clc

do_plot=true;

scale=100;

[db,varlist0]=create_dataset(scale,do_plot);

varlist=fieldnames(varlist0);

%% set up Markov chains

markov_chains=struct('name','syncvol',...
    'number_of_states',3,...
    'controlled_parameters',{{'s'}},...
    'endogenous_probabilities',[],...
    'probability_parameters',[]);

%% Create the VAR

clc

nlags=2;

exog={};

constant=true;

panel=[];

sv0=svar(varlist,exog,nlags,constant,panel,markov_chains);

%% set up restrictions

% syntax is alag(eqtn,vname)
%-------------------------------
lin_restr={
    % first equation or "FFR" equation
    %----------------------------------
    'a1(1,pi)=0'
    'a2(1,pi)=0'
    'a1(1,ygap)=0'
    'a2(1,ygap)=0'
    'a2(1,FFR)=0'
    % second equation or "pi" equation
    %----------------------------------
    'a0(2,FFR)=0'
    'a1(2,FFR)=0'
    'a2(2,FFR)=0'
    'a1(2,ygap)=0'
    'a2(2,ygap)=0'
    % third equation or "ygap" equation
    %-----------------------------------
    'a1(3,FFR)=0'
    'a2(3,FFR)=0'
    'a1(3,pi)=0'
    'a2(3,pi)=0'
    'a0(3,pi)+a0(3,FFR)=0'
    };
nonlin_restr={
    'a0(3,FFR)>=0'
    'a1(1,FFR)>=0'
    'a1(1,FFR)<=1'
    };

restrictions=[lin_restr;nonlin_restr];

%% set priors

% priors for the VAR coefficients
%--------------------------------
var_prior=svar.prior_template();
var_prior.type='sz';
% priors for the syncvol transition probabilities
%------------------------------------------------
switch_prior=struct();
switch_prior.dirichlet_1={0.1,'syncvol_tp_1_2',0.2,'syncvol_tp_1_3',0.2};
switch_prior.dirichlet_2={0.1,'syncvol_tp_2_1',0.2,'syncvol_tp_2_3',0.2};
switch_prior.dirichlet_3={0.1,'syncvol_tp_3_1',0.2,'syncvol_tp_3_2',0.2};

prior=struct();

prior.var=var_prior;

prior.nonvar=switch_prior;

%% Find posterior mode
clc

sv=estimate(sv0,db,{'1960Q1','2015Q2'},prior,restrictions,'fmincon');

%% Find posterior mode : known regimes
clc

db2=db;

nobs=get(db.FFR,'NumberOfObservations');

h=sv0.nregs;

db2.hist_regimes=ts(get(db.FFR,'start'),randi(h,nobs,1));

is_fixed_regime=true;

sv2=estimate(sv0,db2,{'1960Q1','2015Q2'},prior,restrictions,'fmincon',is_fixed_regime);

%% estimates
clc

pmode=posterior_mode(sv)

%% Printing estimates
clc

print_structural_form(sv)

%% Printing solution
clc

print_solution(sv)

%% plot smoothed state and regime probabilities
clc

close all

plot_probabilities(sv)

%% plots probabilities against data
close all

plot_data_against_probabilities(sv,'regime')

%% Impulse responses

myirfs=irf(sv);

snames=fieldnames(myirfs);

figure('name','Simple impulse response functions')

iter=0;

jter=0;

for ishk=1:numel(snames)
    
    shk=snames{ishk};
    
    for iv=1:numel(varlist)
        
        v=varlist{iv};
        
        iter=iter+1;
        
        subplot(3,3,iter)
        
        plot(myirfs.(shk).(v),'linewidth',2)
        
        if iter<=3,title(v),end
        
        if any(iter==[1,4,7])
            
            jter=jter+1;
            
            ylabel([shk,'(',varlist{jter},' shock)'])
            
        end
        
    end
    
end

%% Generalized impulse responses

shock_names=[];

irf_periods=[];

params=[];

girf_setup=struct();% girf_setup=struct('nsims',300);

mygirfs=irf(sv,shock_names,irf_periods,params,girf_setup);

figure('name','Generalized impulse response functions')

iter=0;

jter=0;

for ishk=1:numel(snames)
    
    shk=snames{ishk};
    
    for iv=1:numel(varlist)
        
        v=varlist{iv};
        
        iter=iter+1;
        
        subplot(3,3,iter)
        
        plot(mygirfs.(shk).(v),'linewidth',2)
        
        if iter<=3,title(v),end
        
        if any(iter==[1,4,7])
            
            jter=jter+1;
            
            ylabel([shk,'(',varlist{jter},' shock)'])
            
        end
        
    end
    
end

%% Posterior sampling

[ff,lb,ub,x0,vcov,self]=pull_objective(sv);

options=struct();
options.alpha=0.234;
options.thin=10;
options.burnin=10^3;
options.N=2*10^4;
options.nchain=2;

lb(~isfinite(lb))=-500;
ub(~isfinite(ub))=500;

results=mh_sampler(ff,lb,ub,options,x0,vcov);

%% MCMC diagnostics
clc

pnames=fieldnames(pmode);

a2tilde_to_a=sv.estim_.linres.a2tilde_to_a;

mcmcobj=mcmc(results.pop,pnames,0,1,1,a2tilde_to_a);

myList={'syncvol_tp_1_2','syncvol_tp_1_3','syncvol_tp_2_1','syncvol_tp_2_3',...
    'syncvol_tp_3_1','syncvol_tp_3_2','s_1_1_syncvol_1','s_1_1_syncvol_2',...
    's_1_1_syncvol_3','s_2_2_syncvol_1','s_2_2_syncvol_2','s_2_2_syncvol_3',...
    's_3_3_syncvol_1','s_3_3_syncvol_2','s_3_3_syncvol_3'};

chain_id=[];

for ii=1:numel(myList)
    
    v=myList{ii};
    
    figure('name',['Univariate convergence diagnostics for ',v]);
    
    subplot(3,2,1),autocorrplot(mcmcobj,v,chain_id,40); title('autocorrelogram')
    
    subplot(3,2,2),densplot(mcmcobj,v,chain_id,250); title('density')
    
    subplot(3,2,3),meanplot(mcmcobj,v,chain_id); title('recursive mean')
    
    subplot(3,2,4),psrf_plot(mcmcobj,v); title('pot. scale red. factor')
    
    subplot(3,2,5),traceplot(mcmcobj,v,chain_id,20); title('trace')
    
    pause
    
    close
    
end


%% Marginal data density
clc

mddobj=mdd(results,ff,lb,ub);

fprintf('Importance sampling::%0.2f\n',is(mddobj,[],mdd.global_options))

fprintf('Reciprocal Importance sampling::%0.2f\n',ris(mddobj,[],mdd.global_options))

fprintf('Meng and Wong''s bridge::%0.2f\n',bridge(mddobj,true,mdd.global_options))

fprintf('Ulrich Mueller::%0.2f\n',mueller(mddobj,[],mdd.global_options))

fprintf('Laplace::%0.2f\n',laplace(mddobj))

fprintf('Sims, Waggoner and Zha::%0.2f\n',swz(mddobj,[],mdd.global_options))

fprintf('Laplace MCMC::%0.2f\n',laplace_mcmc(mddobj))

fprintf('Chib and Jeliazkov::%0.2f\n',cj(mddobj,[],mdd.global_options))

%% Out-of sample forecasts at the mode
shock_uncertainty=false;

nsteps=12;

mycast=forecast(sv,[],[],[],nsteps,shock_uncertainty);

do_plot_unconditional_forecasts(mycast,sv)

%% Conditional forecast on ygap

% conditional information
%-------------------------
ygap=scale*(-0.025:0.005:8*0.005).';
cond_db=struct('ygap',ts('2015Q3',ygap));

% options for the exercise
%--------------------------
myoptions=struct('cbands',[10,20,50,80,90],'do_plot',true,'nsteps',20,...
    'param_uncertainty',false,'shock_uncertainty',true,'ndraws',200);

% do it
%-------
tic
[fkst,bands,hdl]=do_conditional_forecasts(mest,db,cond_db,Results.pop,myoptions);
fprintf('\n\n Computing conditional forecasts took %0.4f minutes\n\n',toc/60);
%% Median conditional forecasts
% we take advantage of the fact that in the bands above we specified 50 in
% the bands above 
%--------------------------------------------------------------------------
figure('name','Median Conditional Forecasts');
nvars=mest.endogenous.number;
for ivar=1:nvars
    thisname=mest.endogenous.name{ivar};
    subplot(nvars,1,ivar)
    plot(bands.(thisname)('ci_50'),'linewidth',2)
    title(mest.endogenous.tex_name{ivar})
end


%% irfs distributions

clc

shock_names=[];

params=[results.pop.x];

myirfs2=irf(sv,shock_names,[],params);

%% plotting of fan chargs

for ireg=1:sv.nregs
    
    reg=sprintf('regime_%d',ireg);
    
    figure('name',['Simple impulse response functions(',reg,')'])
    
    iter=0;
    
    jter=0;
    
    for ishk=1:numel(snames)
        
        shk=snames{ishk};
        
        for iv=1:numel(varlist)
            
            v=varlist{iv};
            
            iter=iter+1;
            
            subplot(3,3,iter)
            
            ffr=fanchart(myirfs2.(shk).(v),[30,50,70,90]);
            
            plot_fanchart(ffr.(reg),'r')
            
            if iter<=3,title(v),end
            
            if any(iter==[1,4,7])
                
                jter=jter+1;
                
                ylabel([shk,'(',varlist{jter},' shock)'])
                
            end
            
        end
        
    end
    
end

