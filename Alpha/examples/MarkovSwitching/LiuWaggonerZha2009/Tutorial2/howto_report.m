% housekeeping
close all
clear all
clc

%% initialize the report with default options
myabstract={['This report investigates switches in the parameters of a ',...
    ' simple new Keynesian model estimated on US data. 4 variants of the model',...
    'are estimated: (i) the first variant allows for switches in volalitily only; ',...
    '(ii) the second model allows for switches in the policy parameters only; ',...
    '(iii) the third specification allows for switches in both policy parameters ',...
    'and the volatility of shocks, but in a synchronized fashion; (iv) the ',...
    'third variant allows for independent switches in both parameters and the ',...
    'volatility of shocks.']
    ' ' % leave a blank space to start the next sentence in a new line
    'We find ample evidence in favor of switching parameters... '
    };
mytitlepage=struct('title','Switches in the US Macroeconomic Data: Policy or Volatility?',...
    'address','Your address',...
    'author','Your name',...
    'email','your email',...
    'date',datestr(now),...
    'abstract',{myabstract});
xrep=rise_report('report_name','newkeynesswitch',...
    'titlepage',mytitlepage);
xrep.pagebreak();
%% What to include in the report
add_declaration_description=true;
add_model_equations=true;
add_estimation_results=true;
add_duration_probabilities=true;
add_data_plot=false;
add_data_against_transition_probabilities=true;
add_irfs=true;
add_historical_decomposition=true;
add_counterfactual=false;
add_autocorrelations=true;
add_variance_decomposition=true;
add_forecast=true;
add_shock_correlations=true;
add_empirical_distribution_of_shocks=true;
add_smoothed_shocks=true;
add_shocks_correlations=true;

%% maximum number of rows and cols per figure
graph_nrows=4;
graph_ncols=3;
%% variables of interest
irf_shocks={'ED','ER','ES'};
my_varlist={'PAI','R','X'};
mystates={'coef_2','vol_2'};
mystate_labels={'low monetary policy response regime','High volatility regime'};
%% load the estimated models and relabel them if necessary
nk_models={'volatilityOnly','policyOnly','volatilityPolicySame','volatilityPolicyIndependent'};
nk_models_newlabels={'volOnly','polOnly','volPolSame','volPolInd'};

model_objects=rise.empty(0);
for imod=1:numel(nk_models)
    tmp=load([nk_models{imod},filesep,'estimation',filesep,'estimated_model']);
    obj=tmp.obj;
    model_objects(imod,1)=obj;
    model_objects(imod,1)=set_properties(model_objects(imod,1),'filename',nk_models_newlabels{imod});
end

%% locate the variables of interest
endo_locs=locate_variables(my_varlist,{model_objects(1).varendo.name});
var_list_tex_names={model_objects(1).varendo(endo_locs).tex_name};
graph_nstar=graph_nrows*graph_ncols;
nvar=numel(my_varlist);
n_myvar_figs=ceil(nvar/graph_nstar);
%% get the list of the observed variables
obsList={model_objects(1).varobs.name};
obsListTexnames={model_objects(1).varobs.tex_name};
nvarobs=numel(obsList);
n_myobs_figs=ceil(nvarobs/graph_nstar);
%% exogenous
shock_list={model_objects(1).varexo.name};
shock_list_tex_names={model_objects(1).varexo.tex_name};
nshocks=numel(shock_list);
n_myshocks_figs=ceil(nshocks/graph_nstar);

%% plot the data

if add_data_plot
    db=pages2struct(model_objects(1).options.data);
    add_str='';
    for ifig=1:n_myobs_figs
        if n_myobs_figs>1
            add_str=['(',int2str(ifig),')'];
        end
        fig_title=['Observed data from the US',add_str];
        fig=figure('name',fig_title);
        [Remains,r,c]=number_of_rows_and_columns_in_figure(ifig,nvarobs,graph_nrows,graph_ncols);
        for ivar=(ifig-1)*graph_nstar+1:min(nvarobs,ifig*graph_nstar)
            ii=ivar-(ifig-1)*graph_nstar;
            subplot(r,c,ii)
            plot_window(db.(obsList{ivar}),'',[],@plot,'linewidth',2);
            title([obsListTexnames(ivar),obsList(ivar)],'fontsize',12);
        end
        myfigure=struct('name',fig,'title',fig_title,'scale',0.85);
        figure(xrep,myfigure);
        xrep.pagebreak();
    end
end

%% declarations
if add_declaration_description
    % here we are using the last model has it contains all the parameters,
    % unlike say the first one.
    % endogenous
    %-----------
    mytable=struct('title','Declarations: Legend for endogenous variables',...
        'table',{report(model_objects(end),'rep_type','varendo')});
    table(xrep,mytable);
    % exogenous
    %-----------
    mytable=struct('title','Declarations: Legend for exogenous variables',...
        'table',{report(model_objects(end),'rep_type','varexo')});
    table(xrep,mytable);
    xrep.pagebreak();
    % parameters
    %-----------
    mytable=struct('title','Declarations: Legend for parameters',...
        'table',{report(model_objects(end),'rep_type','parameters')});
    table(xrep,mytable);
    xrep.pagebreak();
    % parameters
    %-----------
    mytable=struct('title','Declarations: Legend for observables',...
        'table',{report(model_objects(end),'rep_type','varobs')});
    table(xrep,mytable);
    xrep.pagebreak();

end
%% add the model equations
if add_model_equations
    myverbatim=struct('title','Model equations',...
        'list',{report(model_objects(end),'rep_type','equations')});
    verbatim(xrep,myverbatim);
    xrep.pagebreak();
end
%% add the estimation results
if add_estimation_results
    mytable=struct('title','Posterior modes',...
        'table',{report(model_objects,'rep_type','estimation')},...
        'longtable',true);
    table(xrep,mytable);
    xrep.pagebreak();
    mytable=struct('title','Estimation statistics',...
        'table',{report(model_objects,'rep_type','estimation_statistics')});
    table(xrep,mytable);
    xrep.pagebreak();
end
%% add duration probabilities

if add_duration_probabilities
for imod=1:numel(nk_models)%nk_models
    mytable={'','Probability of staying','Duration (quarters)'};
    est_par_names={model_objects(imod).estimation.priors.name};
    for ich=1:numel(model_objects(imod).markov_chains)
        chain= model_objects(imod).markov_chains(ich).name;
        mytable=[mytable;{['markov chain name=',chain],'',''}]; %#ok<*AGROW>
        for ist=1:numel(model_objects(imod).markov_chains(ich).states)
            guide=[chain,'_tp_',int2str(ist)];
            lg=length(guide);
            locs= strncmp(guide,est_par_names,lg);
            sumProbs=sum([model_objects(imod).estimated_parameters(locs).mode]);
            duration=1/sumProbs;
            mytable=[mytable;{['state ',int2str(ist)],1-sumProbs,duration}];
        end
    end
    mytable=struct('title',[nk_models{imod},' model: Expected duration in each state'],...
        'table',{mytable});
    table(xrep,mytable);
    xrep.pagebreak();
    
    if model_objects(imod).markov_chains.regimes_number>2
        Q=model_objects(imod).Q;
        qsize=size(Q,1);
        Regimes=model_objects(imod).Regimes;
        mytable={'','Probability of staying','Duration (quarters)'};
        for iq=1:qsize
            thisprobs=Q(iq,:);thisprobs(iq)=[];
            sumProbs=sum(thisprobs);
            duration=1/sumProbs;
            this_regime='';
            for ich=1:numel(model_objects(imod).markov_chains)
                chain= model_objects(imod).markov_chains(ich).name;
                this_regime=[this_regime,', ',chain,'=',int2str(Regimes(iq,ich))];
            end
            this_regime=strtrim(this_regime(2:end));
            mytable=[mytable;{['regime ',int2str(iq),'(',this_regime,')'],1-sumProbs,duration}];
        end
        mytable=struct('title',[nk_models{imod},' model: Expected duration in each Regime'],...
            'table',{mytable});
        table(xrep,mytable);
        xrep.pagebreak();
    end
end
end
%% plot the observables against the transition probabilities

if add_data_against_transition_probabilities
    for imod=1:numel(model_objects)
        if model_objects(imod).markov_chains.regimes_number>1
            thisstates=mystates;
            thislabels=mystate_labels;
            discard=false(1,numel(thisstates));
            for ii=1:numel(thisstates)
                discard(ii)=~ismember(thisstates{ii},model_objects(imod).state_names);
            end
            thisstates=thisstates(~discard);
            thislabels=thislabels(~discard);
            nstates=numel(thisstates);
            %----------------------------
            mytitle=[nk_models{imod},' model: Smoothed probabilities'];
            fig=figure('name',mytitle);
            for istate=1:nstates
                subplot(nstates,1,istate)
                plot(model_objects(imod).Filters.smoothed_probabilities.(thisstates{istate}),...
                    'linewidth',2)
                title([thislabels{istate},'(chain: ',...
                    thisstates{istate}(1:end-2),', state: ',...
                    thisstates{istate}(end),')'],'fontsize',12)
            end
            myfigure=struct('name',fig,'title',mytitle,'angle',90);
            figure(xrep,myfigure);
            xrep.pagebreak();
            %----------------------------
            smoothed=model_objects(imod).Filters.Expected_smoothed_variables;
            for istate=1:nstates
                highvol=model_objects(imod).Filters.smoothed_probabilities.(thisstates{istate});
                add_str='';
                for ifig=1:n_myobs_figs
                    if n_myobs_figs>1
                        add_str=['(',int2str(ifig),')'];
                    end
                    fig_title=[nk_models{imod},' model: Observed series against ',...
                        thislabels{istate},' ',add_str];
                    fig=figure('name',fig_title);
                    [Remains,r,c]=number_of_rows_and_columns_in_figure(ifig,nvarobs,graph_nrows,graph_ncols);
                    for ivar=(ifig-1)*graph_nstar+1:min(nvarobs,ifig*graph_nstar)
                        vname=obsList{ivar};
                        ii=ivar-(ifig-1)*graph_nstar;
                        subplot(r,c,ii)
                        plotyy(smoothed.(vname),highvol,'linewidth',2)
                        title(obsListTexnames{ivar},'fontsize',12)
                    end
                    myfigure=struct('name',fig,'title',fig_title,'scale',0.85);
                    figure(xrep,myfigure);
                    xrep.pagebreak();
                end
            end
        end
    end
end
%% compute impulse responses
if add_irfs
    myirfs=irf(model_objects,'irf_periods',10,'irf_type','girf');
    for ishock=1:numel(shock_list)
        shock=shock_list{ishock};
        if ~ismember(shock,irf_shocks)
            continue
        end
        add_str='';
        for ifig=1:n_myvar_figs
            if n_myvar_figs>1
                add_str=['(',int2str(ifig),')'];
            end
            fig_title=['IRFs to a ',shock_list_tex_names{ishock}, ' shock ',add_str];
            fig=figure('name',fig_title);
            [Remains,r,c]=number_of_rows_and_columns_in_figure(ifig,nvar,graph_nrows,graph_ncols);
            for ivar=(ifig-1)*graph_nstar+1:min(nvar,ifig*graph_nstar)
                endovar=my_varlist{ivar};
                ii=ivar-(ifig-1)*graph_nstar;
                subplot(r,c,ii)
                plot(myirfs.(shock).(endovar),'linewidth',2)
                title(var_list_tex_names{ivar},...
                    'fontsize',12)
                if ivar==1
                    leg=legend(nk_models);
                    set(leg,'interpreter','none')
                end
            end
            orient(fig,'tall')
            myfigure=struct('name',fig,'title',fig_title,'scale',0.85);
            figure(xrep,myfigure);
            xrep.pagebreak();
        end
    end
end
%% historical decomposition of shocks
if add_historical_decomposition
    histdec=historical_decomposition(model_objects);
    for ivar=1:nvar
        if ivar==1
            contrib_names=histdec{1}.(vname).varnames;
            locs=locate_variables(contrib_names,shock_list,true);
            for ii=1:numel(locs)
                if isnan(locs(ii))
                    continue
                end
                contrib_names{ii}=shock_list_tex_names{locs(ii)};
            end
        end
        vname=my_varlist{ivar};
        fig_title=['historical decomposition of ',var_list_tex_names{ivar}];
        fig=figure('name',fig_title);
        for imod=1:numel(nk_models)
            subplot(numel(nk_models),1,imod)
            plot_decomp(histdec{imod}.(vname))
            title(nk_models{imod},'interpreter','none',...
                'fontsize',12)
        end
        orient(fig,'tall')
        hleg=legend(contrib_names,'location','SW','orientation','horizontal');
        myfigure=struct('name',fig,'title',fig_title,'scale',0.85);
        figure(xrep,myfigure);
        xrep.pagebreak();
    end
end

%% counterfactual: what if the economy had been ...
if add_counterfactual
    for ii=1:10
        disp('the code for doing counterfactuals has been temporarily removed')
    end
end
%% variance decomposition

if add_variance_decomposition
    vardec=variance_decomposition(model_objects);
    % plot the decomposition
    for ivar=1:nvarobs
        vname=obsList{ivar};
        fig_title=['Variance decomposition of ',obsListTexnames{ivar}];
        fig=figure('name',fig_title);
        for imod=1:numel(nk_models)
            h0=subplot(numel(nk_models),1,imod);
            tmp=100*vardec{imod}.conditional.(vname);
            tmp.varnames=vardec{imod}.conditional.(vname).varnames;
            plt=plot_window(tmp,1,24,@plot_decomp);
            title([obsListTexnames{ivar},'(',nk_models{imod},')'],...
                'fontsize',12,'interp','none')
            if imod==1
                contrib_names=vardec{imod}.conditional.(vname).varnames;
                locs=locate_variables(contrib_names,shock_list,true);
                for ii=1:numel(locs)
                    if isnan(locs(ii))
                        continue
                    end
                    contrib_names{ii}=shock_list_tex_names{locs(ii)};
                end
                legend(contrib_names,'location','NorthEast');
            end
        end
        orient(fig,'tall')
        myfigure=struct('name',fig,'title',fig_title,'scale',0.85);
        figure(xrep,myfigure);
        xrep.pagebreak();
    end
end
%% real-time forecasting performance of the estimated model
if add_forecast
    map=getappdata(0,'rise_default_plot_colors');
    
    [ts_fkst,ts_rmse,rmse,Updates]=forecast_real_time(model_objects);
    TimeInfo=ts_fkst{1}.(obsList{1}).TimeInfo;
    pp=[];
    Q=[];
    for iobs=1:nvarobs
        fig_title=['real-time forecasts for ',obsListTexnames{iobs}];
        fig=figure('name',fig_title);
        for imod=1:numel(nk_models)
            subplot(numel(nk_models)+1,1,imod)
            [h,pp]=plot_real_time(ts_fkst{imod}.(obsList{iobs}),pp);
            % add the steady state
            vloc=locate_variables(obsList{iobs},{model_objects(imod).varendo.name});
            hold on
            xlim=get(h,'xlim');
            ergodic_mean=model_objects(imod).varendo(vloc).det_steady_state;
            if numel(ergodic_mean)>1
                Q=model_objects(imod).Q;
                nreg=size(Q,1);
                pai=[eye(nreg)-Q';ones(1,nreg)]\[zeros(nreg,1);1];
                ergodic_mean=ergodic_mean*pai;
            end
            plot(xlim,ergodic_mean*ones(1,2),'color',[1,0,0])
            hold off
            title([obsListTexnames{iobs},'(',nk_models{imod},')'],...
                'fontsize',12,'interpreter','none')
            grid on
        end
        % plot the rmses
        subplot(numel(nk_models)+1,1,numel(nk_models)+1)
        plot(ts_rmse.(obsList{iobs}))
        title('RMSE','fontsize',12)
        leg=legend(nk_models);
        set(leg,'interpreter','none')
        orient(fig,'tall')
        myfigure=struct('name',fig,'title',fig_title,'scale',0.85);
        figure(xrep,myfigure);
        xrep.pagebreak();
    end
end
%% vector autocorrelations
if add_autocorrelations
    [Auto,retcode]=theoretical_autocorrelations(model_objects,'ar',40);
    for ii=1:numel(Auto)
        Auto{ii}=Auto{ii}(endo_locs,endo_locs,:);
    end
    ar=size(Auto{1},3);
    xx=(1:ar)';
    add_str='';
    nfigs=ceil(nvar^2/graph_nstar);
    [first_index,second_index]=ind2sub([nvar,nvar],1:nvar^2);
    for ifig=1:nfigs
        if nfigs>1
            add_str=['(',int2str(ifig),')'];
        end
        fig_title=['Vector autocorrelations',add_str];
        fig=figure('name',fig_title);
        [Remains,r,c]=number_of_rows_and_columns_in_figure(ifig,nvar^2,graph_nrows,graph_ncols);
        for ivar=(ifig-1)*graph_nstar+1:min(nvar^2,ifig*graph_nstar)
            vname1=var_list_tex_names{second_index(ivar)};
            vname2=var_list_tex_names{first_index(ivar)};
            ii=ivar-(ifig-1)*graph_nstar;
            ar_data=[];
            for imod=1:numel(nk_models)
                ar_data=[ar_data,squeeze(Auto{imod}(second_index(ivar),first_index(ivar),:))];
            end
            subplot(r,c,ii)
            plot(xx,ar_data,'linewidth',2)
            axis tight
            if ii==1||rem(ii,c)==1
                ylabel(vname1)
            end
            if ii>(r-1)*c
                xlabel(vname2)
            end
            if ii==1
                legend(nk_models)
            end
        end
        orient(fig,'tall')
        myfigure=struct('name',fig,'title',fig_title,'scale',0.85);
        figure(xrep,myfigure);
        xrep.pagebreak();
    end
end
%% empirical distribution of shocks

if add_smoothed_shocks
    add_str='';
    for ifig=1:n_myshocks_figs
        if n_myshocks_figs>1
            add_str=['(',int2str(ifig),')'];
        end
        fig_title=['Smoothed shocks',add_str];
        fig=figure('name',fig_title);
        [Remains,r,c]=number_of_rows_and_columns_in_figure(ifig,nshocks,graph_nrows,graph_ncols);
        for ivar=(ifig-1)*graph_nstar+1:min(nshocks,ifig*graph_nstar)
            ii=ivar-(ifig-1)*graph_nstar;
            vshock=shock_list{ii};
            myshock_data=model_objects(1).Filters.Expected_smoothed_shocks.(vshock);
            for imod=2:numel(model_objects)
                myshock_data=[myshock_data,model_objects(imod).Filters.Expected_smoothed_shocks.(vshock)];
            end
            subplot(r,c,ii)
            plot(myshock_data,'linewidth',2)
            axis tight
            title(shock_list_tex_names{ivar},...
                'fontsize',12)
            if ii==1
                leg=legend({model_objects.filename});
                set(leg,'interpreter','none')
            end
        end
        orient(fig,'tall')
        myfigure=struct('name',fig,'title',fig_title,'scale',0.85);
        figure(xrep,myfigure);
        xrep.pagebreak();
    end
end

%% correlation of shocks
if add_shocks_correlations
    for imod=1:numel(model_objects)
        shock_corr=corrcoef(rise_time_series.collect(model_objects(imod).Filters.Expected_smoothed_shocks));
        mytable=[
            [' ',shock_list]
            shock_list',num2cell(shock_corr)
            ];
        mytable=struct('title',[nk_models{imod},' model: Shock correlation structure '],...
            'table',{mytable});
        table(xrep,mytable);
        xrep.pagebreak();
    end
end
%% empirical distribution of shocks

if add_empirical_distribution_of_shocks
    add_str='';
    for ifig=1:n_myshocks_figs
        if n_myshocks_figs>1
            add_str=['(',int2str(ifig),')'];
        end
        fig_title=['Empirical PDF of shocks',add_str];
        fig=figure('name',fig_title);
        [Remains,r,c]=number_of_rows_and_columns_in_figure(ifig,nshocks,graph_nrows,graph_ncols);
        for ivar=(ifig-1)*graph_nstar+1:min(nshocks,ifig*graph_nstar)
            ii=ivar-(ifig-1)*graph_nstar;
            vshock=shock_list{ii};
            myshock_data=cell(1,numel(model_objects));
            lb=inf;
            ub=-inf;
            for imod=1:numel(model_objects)
                myshock_data{imod}=double(model_objects(imod).Filters.Expected_smoothed_shocks.(vshock));
                lb=min(lb,min(myshock_data{imod}));
                ub=max(ub,max(myshock_data{imod}));
            end
            myff=[];
            for imod=1:numel(model_objects)
                [ff_,xx_]=distributions.kernel_density(myshock_data{imod},lb,ub,'normal',250);
                myff=[myff,ff_(:)];
            end
            subplot(r,c,ii)
            plot(xx_(:),myff,'linewidth',2)
            axis tight
            title(shock_list_tex_names{ivar},'fontsize',12)
            if ii==min(nshocks,ifig*graph_nstar)
                leg=legend(nk_models,'location','SO','orientation','horizontal');
                set(leg,'interpreter','none')
            end
        end
        orient(fig,'tall')
        myfigure=struct('name',fig,'title',fig_title,'scale',0.85);
        figure(xrep,myfigure);
        xrep.pagebreak();
    end
end

%% now create the report
publish(xrep)
close all

