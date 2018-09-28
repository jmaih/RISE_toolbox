%% housekeeping
clear
close all
clc
%% RISE the model
m=rise('nkp');

%% calibrated parameters and priors

start_from_mode=false;

[p,priors]=create_parameters(start_from_mode);

m=set(m,'parameters',p);

%% data

[data]=create_data();

%% estimate model
clc
log_vars=[];%{'A','C','G','LAMBDA','PAI','Q','R','THETA','X','Y','Z'};
[mest,filtration]=estimate(m,'data',data,'estim_priors',priors,...
    'solve_log_approx_vars',log_vars);
%% Impulse response functions
myirfs=irf(mest);
% % close,shock='EPS_Z';plot(cumsum(myirfs.(shock).Z)+myirfs.(shock).Y)
close all
% mest=set(mest,'tex_name',{
%     'Y','Output'
%     'PAI','Inflation'
%     'X','Output gap'
%     'R','Interest rate'
%     'EPS_E','Cost push'
%     'EPS_R','monetary policy'
%     'EPS_Z','Technology'
%     'EPS_A','Preference'
%     });
tex=get(m,'tex');

figure('name','Impulse responses');
myvars={'Y','PAI','R','X'};
nvars=numel(myvars);
iter=0;
for ivar=1:nvars
    vname=myvars{ivar};
    vtex=tex.(vname);
    for ishock=1:mest.exogenous.number(1)
        shock=mest.exogenous.name{ishock};
        shock_tex=tex.(shock);
        iter=iter+1;
        subplot(4,4,iter)
        plot('0:20',100*myirfs.(shock).(vname),'linewidth',2)
        title([vtex,' to ',shock_tex])
    end
end
orient landscape
%% variance decomposition
clc
mydec=variance_decomposition(mest);
myvars={'G','PAI','R','X'};
for ivar=1:numel(myvars)
    vname=myvars{ivar};
    vtex=tex.(vname);
    dec=mydec.conditional.(vname)('1,4,8,12,20,40');
    if ivar==1
        shock_pos=locate_variables(mest.exogenous.name,dec.varnames);
    end
    dec=100*dec;
    dec.varnames=mest.exogenous.tex_name(shock_pos);
    fprintf('\n\n %s \n', vtex);
    display(dec)
end
%% Shocks analysis
clc
recessions={'1990Q3:1992Q3','2001Q1:2003Q1','2007Q4:2009Q4'};
pest=get(mest,'parameters');
nshocks=mest.exogenous.number(1);
sigmas=cell(1,nshocks);
for irecess=1:numel(recessions)
    recess=recessions{irecess};
    db=ts.empty(0);
    for ishock=1:nshocks
        shock=mest.exogenous.name{ishock};
        last_letter=lower(shock(end));
        sig=eval(['pest.sig_',last_letter]);
        db=[db,sig*filtration.smoothed_shocks.(shock)(recess)];
        if irecess==1
            sigmas{ishock}=sig;
        end
    end
    db.varnames=mest.exogenous.tex_name;
    fprintf('\n\n %s recession \n', recess);
    display(db)
    if irecess==numel(recessions)
        fprintf('\n\n Estimated standard deviations \n');
        disp([db.varnames;sigmas])
    end
end
%% Counterfactuals
clc
close all
db=filtration.smoothed_variables;
dbshocks=filtration.smoothed_shocks;
db=pages2struct(ts.collect(db,dbshocks));
for irecess=1:numel(recessions)
    recess=recessions{irecess};
    serials=date2serial(recess);
    nsteps=numel(serials)-1;
    figure('name',['Counterfactual output paths for ',recess,' recession'])
    for ishock=1:nshocks
        shock=mest.exogenous.name{ishock};
        fkst=forecast(mest,'data',pageify(serials(2)-1,db),...
            'forecast_start_date',serials(2),...
            'forecast_cond_exo_vars',shock,...
            'forecast_nsteps',nsteps);
        subplot(2,2,ishock)
        toplot=[db.GHAT+log(pest.zss),fkst.GHAT+log(pest.zss)];
        toplot=cumsum(toplot(recess));
        benchmark=double(toplot);
        benchmark=max(benchmark(:));
        plot(toplot,'linewidth',2)%plot((toplot-benchmark)/benchmark,'linewidth',2)
        title(mest.exogenous.tex_name{ishock})
        if ishock==1
            legend({'actual','forecast'})
        end
    end
end
%% full sample estimates of monetary policy shocks
figure('name','Full-sample estimates of monetary policy shocks')
plot(filtration.smoothed_shocks.EPS_R,'linewidth',2)