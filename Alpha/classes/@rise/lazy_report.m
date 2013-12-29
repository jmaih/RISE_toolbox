function out=lazy_report(obj,varargin)

if isempty(obj)
    out=struct('lazy_report_abstract','',...
        'lazy_report_title','untitled',...
        'lazy_report_authors','',...
        'lazy_report_name','',...
        'lazy_report_date','',...
        'lazy_report_shock_list','',...
        'lazy_report_var_list','');
    return
end

obj=set_options(obj,varargin{:});

out=[];
%% maximum number of rows and cols per figure
nmod=numel(obj);
graph_nrows=4;
graph_ncols=3;
start=[];
finish=[];
modelsLegs={obj.legend};
for imod=1:nmod
    if isempty(obj(imod).legend)
        modelsLegs{imod}=sprintf('model #%0.0f',imod);
    end
end
%% variables of interest
% irf_shocks={'ED','ER','ES'};
shock_list=obj(1).options.lazy_report_shock_list;
if isempty(shock_list)
    shock_list=obj(1).exogenous.name;
end
var_list=obj(1).options.lazy_report_var_list;
if isempty(var_list)
    var_list=obj(1).endogenous.name(obj(1).endogenous.is_original & ~obj(1).endogenous.is_auxiliary);
end
% mystates={'coef_2','vol_2'};
% mystate_labels={'low monetary policy response regime','High volatility regime'};


%% What to include in the report
add_declaration_description=true;
add_model_equations=true;
add_estimation_results=true;
add_duration_probabilities=true;
add_data_plot=true;
add_transition_probabilities=true;
add_data_against_transition_probabilities=true;
add_irfs=true;
add_historical_decomposition=true;
add_counterfactual=false;
add_autocorrelations=true;
add_variance_decomposition=true;
add_forecast=true;
add_empirical_distribution_of_shocks=true;
add_smoothed_shocks=true;
add_shocks_correlations=true;
add_estimated_shocks=true;

%% locate the variables of interest
endo_locs=locate_variables(var_list,obj(1).endogenous.name);
var_list_tex_names=obj(1).endogenous.tex_name(endo_locs);
graph_nstar=graph_nrows*graph_ncols;
nvar=numel(var_list);
n_myvar_figs=ceil(nvar/graph_nstar);
%% get the list of the observed variables
is_endogenous=obj(1).observables.is_endogenous;
obsList=obj(1).observables.name(is_endogenous);
obsListTexnames=obj(1).observables.tex_name(is_endogenous);
nvarobs=obj(1).observables.number(1);
n_myobs_figs=ceil(nvarobs/graph_nstar);
%% exogenous
exo_locs=locate_variables(shock_list,obj(1).exogenous.name);
shock_list_tex_names=obj(1).exogenous.tex_name(exo_locs);
nshocks=numel(shock_list);
n_myshocks_figs=ceil(nshocks/graph_nstar);
%% collect the smoothed series
smoothed=[];
if ~isempty(obj(1).filtering)
smoothed=cell(1,nmod);
for imod=1:nmod
    smoothed{imod}=mergestructures(obj(imod).filtering.Expected_smoothed_variables,...
    obj(imod).filtering.Expected_smoothed_shocks,...
    obj(imod).filtering.smoothed_regime_probabilities,...
    obj(imod).filtering.smoothed_state_probabilities);
end
end
%% initialize the report 
mytitlepage=struct();
if ~isempty(obj(1).options.lazy_report_title)
    mytitlepage.title=char(obj(1).options.lazy_report_title);
end
    authors=obj(1).options.lazy_report_authors;
if ~isempty(authors)
    if ischar(authors)
        ampersand=find(authors =='&');
        if ~isempty(ampersand)
            myauthors=cell(1,numel(ampersand)+1);
            myauthors{1}=authors(1:ampersand(1)-1);
            for ia=1:numel(ampersand)-1
                myauthors{ia+1}=authors(ampersand(ia)+1:ampersand(ia+1)-1);
            end
            myauthors{end}=authors(ampersand(end)+1:end);
        end
        authors=myauthors;
    end
    mytitlepage.author=authors;
end
lazy_report_date=obj(1).options.lazy_report_date;
if isempty(lazy_report_date)
    lazy_report_date=datestr(now);
end
mytitlepage.date=lazy_report_date;
abstract=obj(1).options.lazy_report_abstract;
if ~isempty(abstract)
    mytitlepage.abstract=char(abstract);
end
lazy_report_name=obj(1).options.lazy_report_name;
if isempty(lazy_report_name)
    dd=dir();
    dd=regexp({dd.name},'rise_report_\d+\.pdf','match');
    dd=[dd{:}];
    if isempty(dd)
        lazy_report_name='rise_report_1';
    else
        highest=0;
        for ifile=1:numel(dd)
            index=str2double(dd{ifile}(13:end));
            highest=max(highest,index);
        end
        lazy_report_name=['rise_report_',sprintf('%0.0f',highest)];
    end
end
xrep=rise_report('report_name',lazy_report_name,'titlepage',mytitlepage);
xrep.pagebreak();
%% plot the data

if add_data_plot
    db=pages2struct(obj(1).options.data);
    add_str='';
    for ifig=1:n_myobs_figs
        if n_myobs_figs>1
            add_str=['(',int2str(ifig),')'];
        end
        fig_title=['Observed data from the US',add_str];
        fig=figure('name',fig_title);
        [~,r,c]=number_of_rows_and_columns_in_figure(ifig,nvarobs,graph_nrows,graph_ncols);
        for ivar=(ifig-1)*graph_nstar+1:min(nvarobs,ifig*graph_nstar)
            ii=ivar-(ifig-1)*graph_nstar;
            subplot(r,c,ii)
            plot_window(db.(obsList{ivar}),'',[],@plot,'linewidth',2);
            title([obsListTexnames(ivar),obsList(ivar)],'fontsize',12);
        end
        xrotate(90)
        myfigure=struct('name',fig,'title',fig_title,'scale',0.85);
        figure(xrep,myfigure);
        xrep.pagebreak();
    end
end

%% declarations
if add_declaration_description
    % endogenous
    %-----------
    mytable=struct('title','Declarations: Legend for endogenous variables',...
        'table',{report(obj(1),'rep_type','endogenous')});
    table(xrep,mytable);
    % exogenous
    %-----------
    mytable=struct('title','Declarations: Legend for exogenous variables',...
        'table',{report(obj(1),'rep_type','exogenous')});
    table(xrep,mytable);
    xrep.pagebreak();
    % parameters
    %-----------
    mytable=struct('title','Declarations: Legend for parameters',...
        'table',{report(obj(1),'rep_type','parameters')});
    table(xrep,mytable);
    xrep.pagebreak();
    % parameters
    %-----------
    mytable=struct('title','Declarations: Legend for observables',...
        'table',{report(obj(1),'rep_type','observables')});
    table(xrep,mytable);
    xrep.pagebreak();

end
%% - the model code and declarations
%------------------------------------
if add_model_equations
    xrep.section('Model code',latex_model_file(obj(1),true,true))
    xrep.pagebreak();
end
%% add the estimation results
%------------------------------------
if add_estimation_results && ~isempty(obj.estimation)
    mytable=struct('title','Posterior modes',...
        'table',{report(obj,'rep_type','estimation')},...
        'longtable',true);
    table(xrep,mytable);
    xrep.pagebreak();
    mytable=struct('title','Estimation statistics',...
        'table',{report(obj,'rep_type','estimation_statistics')});
    table(xrep,mytable);
    xrep.pagebreak();
end

%%
if add_transition_probabilities && ~isempty(obj.filtering)
    for imod=1:nmod
        model=obj(imod);
        % 3-a  The smoothed regime probabilities
        %---------------------------------------------------
        packages=model.markov_chains.regime_names;
        fig_title=['Smoothed Regime Probabilities ',model.legend];
        plot_packages(smoothed{imod},graph_nrows,graph_ncols,packages,start,finish,fig_title)
        xrotate(90)
        myfigure=struct('name',gcf,'title',fig_title,'scale',0.85);
        figure(xrep,myfigure);
        xrep.pagebreak();
        
        % 3-b  The smoothed regime probabilities(all in one)
        %---------------------------------------------------
        packages={model.markov_chains.regime_names};
        fig_title=['Smoothed Regime Probabilities in ',model.legend];
        plot_packages(smoothed{imod},graph_nrows,graph_ncols,packages,start,finish,fig_title)
        xrotate(90)
        myfigure=struct('name',gcf,'title',fig_title,'scale',0.85);
        figure(xrep,myfigure);
        xrep.pagebreak();
        
        % 3-c smoothed state probabilities
        %----------------------------------
        packages=model.markov_chains.state_names;
        fig_title=['Smoothed state Probabilities ',model.legend];
        plot_packages(smoothed{imod},graph_nrows,graph_ncols,packages,start,finish,fig_title)
        xrotate(90)
        myfigure=struct('name',gcf,'title',fig_title,'scale',0.85);
        figure(xrep,myfigure);
        xrep.pagebreak();
    end
end

%% add duration probabilities

if add_duration_probabilities && ~isempty(obj.filtering)
for imod=1:nmod
    mytable={'','Probability of staying','Duration (quarters)'};
    params=get(obj(imod),'parameters');
    pnames=fieldnames(params);
    markov_chains=obj(imod).markov_chains;
    regimes=markov_chains.regimes;
    for ich=1:markov_chains.chains_number
        chain= markov_chains.chain_names{ich};
        if strcmp(chain,'const')
            continue
        end
        nstates=max(cell2mat(regimes(2:end,ich+1)));
        mytable=[mytable;{['markov chain name=',chain],'',''}]; %#ok<*AGROW>
        for ist=1:nstates
            guide=[chain,'_tp_',int2str(ist)];
            lg=length(guide);
            locs= find(strncmp(guide,pnames,lg));
            sumProbs=0;
            for iprob=1:numel(locs)
                sumProbs=sumProbs+params.(pnames{locs(iprob)})(1);
            end
            duration=1/sumProbs;
            mytable=[mytable;{['state ',int2str(ist)],1-sumProbs,duration}];
        end
    end
    mytable=struct('title',[obj(imod).legend,' model: Expected duration in each state'],...
        'table',{mytable});
    table(xrep,mytable);
    xrep.pagebreak();
    
    if obj(imod).markov_chains.regimes_number>2
        Q=obj(imod).solution.Q;
        qsize=size(Q,1);
        regimes=obj(imod).markov_chains.regimes;
        mytable={'','Probability of staying','Duration (quarters)'};
        for iq=1:qsize
            thisprobs=Q(iq,:);thisprobs(iq)=[];
            sumProbs=sum(thisprobs);
            duration=1/sumProbs;
            this_regime='';
            for ich=1:markov_chains.chains_number
                chain= markov_chains.chain_names{ich};
                % if strcmp(chain,'const'),continue,end
                this_regime=[this_regime,', ',chain,'=',int2str(regimes{iq+1,ich+1})];
            end
            this_regime=strtrim(this_regime(2:end));
            mytable=[mytable;{['regime ',int2str(iq),'(',this_regime,')'],1-sumProbs,duration}];
        end
        mytable=struct('title',[obj(imod).legend,' model: Expected duration in each Regime'],...
            'table',{mytable});
        table(xrep,mytable);
        xrep.pagebreak();
    end
end
end
%% plot the observables against the transition probabilities

if add_data_against_transition_probabilities && ~isempty(obj.filtering)
    for imod=1:numel(obj)
        mystates=obj(imod).markov_chains.state_names;
        mystate_labels=obj(imod).markov_chains.state_tex_names;
        if obj(imod).markov_chains.regimes_number>1
            thisstates=mystates;
            thislabels=mystate_labels;
            discard=false(1,numel(thisstates));
            for ii=1:numel(thisstates)
                discard(ii)=~ismember(thisstates{ii},obj(imod).markov_chains.state_names);
            end
            thisstates=thisstates(~discard);
            thislabels=thislabels(~discard);
            nstates=numel(thisstates);
            %----------------------------
            mytitle=[obj(imod).legend,' model: Smoothed probabilities'];
            fig=figure('name',mytitle);
            for istate=1:nstates
                subplot(nstates,1,istate)
                plot(smoothed{imod}.(thisstates{istate}),...
                    'linewidth',2)
                title([thislabels{istate},'(chain: ',...
                    thisstates{istate}(1:end-2),', state: ',...
                    thisstates{istate}(end),')'],'fontsize',12)
            end
            myfigure=struct('name',fig,'title',mytitle,'angle',90);
            figure(xrep,myfigure);
            xrep.pagebreak();
            %----------------------------
            for istate=1:nstates
                highvol=smoothed{imod}.(thisstates{istate});
                add_str='';
                for ifig=1:n_myobs_figs
                    if n_myobs_figs>1
                        add_str=['(',int2str(ifig),')'];
                    end
                    fig_title=[obj(imod).legend,' model: Observed series against ',...
                        thislabels{istate},' ',add_str];
                    fig=figure('name',fig_title);
                    [~,r,c]=number_of_rows_and_columns_in_figure(ifig,nvarobs,graph_nrows,graph_ncols);
                    for ivar=(ifig-1)*graph_nstar+1:min(nvarobs,ifig*graph_nstar)
                        vname=obsList{ivar};
                        ii=ivar-(ifig-1)*graph_nstar;
                        subplot(r,c,ii)
                        plotyy(smoothed{imod}.(vname),highvol,'linewidth',2)
                        title(obsListTexnames{ivar},'fontsize',12)
                    end
                    xrotate(90)
                    myfigure=struct('name',fig,'title',fig_title,'scale',0.85);
                    figure(xrep,myfigure);
                    xrep.pagebreak();
                end
            end
        end
    end
end

%% 12- estimated shocks
%--------------------
if add_estimated_shocks && ~isempty(obj.filtering)
    for imod=1:nmod
        model=obj(imod);
        packages=model.exogenous.name;
        fig_title=['Estimated Shocks(All sample) in ',model.legend];
        plot_packages(smoothed{imod},graph_nrows,graph_ncols,packages,start,finish,fig_title,0)
        xrotate(90)
        myfigure=struct('name',gcf,'title',fig_title,'scale',0.85);
        figure(xrep,myfigure);
        xrep.pagebreak();
        
        % 13- Histograms for the estimated shocks
        %-----------------------------------------
        tight=true;
        packages=model.exogenous.name;
        fig_title=['Histograms of shocks in ',model.legend];
        plot_packages(smoothed{imod},graph_nrows,graph_ncols,packages,start,finish,fig_title,[],@hist,tight)
        myfigure=struct('name',gcf,'title',fig_title,'scale',0.85);
        figure(xrep,myfigure);
        xrep.pagebreak();
        
        % 13-b) empirical distribution of shocks
    end
end
%% - Sample correlation of estimated shocks
%-------------------------------------------
if add_shocks_correlations && ~isempty(obj.filtering)
    for imod=1:numel(model)
        shock_corr=corrcoef(rise_time_series.collect(model(imod).filtering.Expected_smoothed_shocks));
        mytable=[
            [' ',shock_list]
            shock_list',num2cell(shock_corr)
            ];
        mytable=struct('title',[model(imod).legend,' model: Shock correlation structure '],...
            'table',{mytable});
        table(xrep,mytable);
        xrep.pagebreak();
    end
end

%% compute impulse responses
if add_irfs
    myirfs=irf(obj);
    for ishock=1:numel(shock_list)
        shock=shock_list{ishock};
        if ~ismember(shock,shock_list)
            continue
        end
        add_str='';
        for ifig=1:n_myvar_figs
            if n_myvar_figs>1
                add_str=['(',int2str(ifig),')'];
            end
            fig_title=['IRFs to a ',shock_list_tex_names{ishock}, ' shock ',add_str];
            fig=figure('name',fig_title);
            [~,r,c]=number_of_rows_and_columns_in_figure(ifig,nvar,graph_nrows,graph_ncols);
            for ivar=(ifig-1)*graph_nstar+1:min(nvar,ifig*graph_nstar)
                endovar=var_list{ivar};
                ii=ivar-(ifig-1)*graph_nstar;
                subplot(r,c,ii)
                plot(myirfs.(shock).(endovar),'linewidth',2)
                title(var_list_tex_names{ivar},...
                    'fontsize',12)
                if ivar==1
                    leg=legend(myirfs.(shock).(endovar).varnames);
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
if add_historical_decomposition && ~isempty(obj.filtering)
    histdec=historical_decomposition(obj);
    cell_flag=iscell(histdec);
    for ivar=1:nvar
        vname=var_list{ivar};
        if ivar==1
            if cell_flag
            contrib_names=histdec{1}.(vname).varnames;
            else
            contrib_names=histdec(1).(vname).varnames;
            end
            locs=locate_variables(contrib_names,shock_list,true);
            for ii=1:numel(locs)
                if isnan(locs(ii))
                    continue
                end
                contrib_names{ii}=shock_list_tex_names{locs(ii)};
            end
        end
        fig_title=['historical decomposition of ',var_list_tex_names{ivar}];
        fig=figure('name',fig_title);
        for imod=1:nmod
            subplot(nmod,1,imod)
            if cell_flag
                plot_decomp(histdec{imod}.(vname))
            else
                plot_decomp(histdec(imod).(vname))
            end
            title(obj(imod).legend,'interpreter','none',...
                'fontsize',12)
        end
        xrotate(90)
        orient(fig,'tall')
        legend(contrib_names,'location','SW','orientation','horizontal');
        myfigure=struct('name',fig,'title',fig_title,'scale',0.85);
        figure(xrep,myfigure);
        xrep.pagebreak();
    end
end

%% counterfactual: what if the economy had been ...
if add_counterfactual && ~isempty(obj.filtering)
    for ii=1:10
        disp('the code for doing counterfactuals has been temporarily removed')
    end
end
%% variance decomposition

if add_variance_decomposition
    vardec=variance_decomposition(obj);
    % plot the decomposition
    for ivar=1:nvarobs
        vname=obsList{ivar};
        fig_title=['Variance decomposition of ',obsListTexnames{ivar}];
        fig=figure('name',fig_title);
        for imod=1:nmod
            subplot(nmod,1,imod)
            if iscell(vardec)
                tmp=100*vardec{imod}.conditional.(vname);
                tmp.varnames=vardec{imod}.conditional.(vname).varnames;
            else
                tmp=100*vardec(imod).conditional.(vname);
                tmp.varnames=vardec(imod).conditional.(vname).varnames;
            end
            plot_window(tmp,1,24,@plot_decomp);
            thismodel=obj(imod).legend;
            if isempty(thismodel)
                thismodel=obj.filename;
            end
            title([obsListTexnames{ivar},'(',thismodel,')'],...
                'fontsize',12,'interp','none')
            if imod==1
                if iscell(vardec)
                contrib_names=vardec{imod}.conditional.(vname).varnames;
                else
                contrib_names=vardec(imod).conditional.(vname).varnames;
                end
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
if add_forecast && ~isempty(obj.filtering)
   
    [ts_fkst,ts_rmse]=forecast_real_time(obj);
    many_models=iscell(ts_fkst);
    pp=[];
    ergodic_mean=cell(1,nmod);
    get(obj(imod),'sstate')
    for iobs=1:nvarobs
        fig_title=['real-time forecasts for ',obsListTexnames{iobs}];
        fig=figure('name',fig_title);
        for imod=1:nmod
            subplot(nmod+1,1,imod)
            if many_models
            [h,pp]=plot_real_time(ts_fkst{imod}.(obsList{iobs}),pp);
            else
            [h,pp]=plot_real_time(ts_fkst.(obsList{iobs}),pp);
            end
            % add the steady state
            vloc=locate_variables(obsList{iobs},obj(imod).endogenous.name);
            hold on
            xlim=get(h,'xlim');
            if iobs==1
                ergodic_mean{imod}=cell2mat(struct2cell(get(obj(imod),'sstate')));
                if size(ergodic_mean{imod},2)>1
                    Q=obj(imod).solution.Q;
                    nreg=size(Q,1);
                    pai=[eye(nreg)-Q';ones(1,nreg)]\[zeros(nreg,1);1];
                    ergodic_mean{imod}=ergodic_mean{imod}*pai;
                end
            end
            plot(xlim,ergodic_mean{imod}(vloc)*ones(1,2),'color',[1,0,0])
            hold off
            thisleg='';
            if nmod>1
                thisleg=['(',modelsLegs{imod},')'];
            end
            title([obsListTexnames{iobs},thisleg],...
                'fontsize',12,'interpreter','none')
            grid on
        end
        % plot the rmses
        subplot(nmod+1,1,nmod+1)
        plot(ts_rmse.(obsList{iobs}))
        title('RMSE','fontsize',12)
        if nmod>1
            leg=legend(modelsLegs);
            set(leg,'interpreter','none')
        end
        orient(fig,'tall')
        myfigure=struct('name',fig,'title',fig_title,'scale',0.85);
        figure(xrep,myfigure);
        xrep.pagebreak();
    end
end
%% vector autocorrelations
if add_autocorrelations
    [Auto,retcode]=theoretical_autocorrelations(obj,'autocorr_ar',40);
    if retcode
        warning('problem detected in the computation of theoretical autocorrelations')
    end
    cellmode=iscell(Auto);
    if cellmode
    for ii=1:numel(Auto)
        Auto{ii}=Auto{ii}(endo_locs,endo_locs,:);
    end
    ar=size(Auto{1},3);
    else
        Auto=Auto(endo_locs,endo_locs,:);
    ar=size(Auto,3);
    end
    
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
        [~,r,c]=number_of_rows_and_columns_in_figure(ifig,nvar^2,graph_nrows,graph_ncols);
        for ivar=(ifig-1)*graph_nstar+1:min(nvar^2,ifig*graph_nstar)
            vname1=var_list_tex_names{second_index(ivar)};
            vname2=var_list_tex_names{first_index(ivar)};
            ii=ivar-(ifig-1)*graph_nstar;
            ar_data=[];
            for imod=1:nmod
                if cellmode
                    this_auto=Auto{imod}(second_index(ivar),first_index(ivar),:);
                else
                    this_auto=Auto(second_index(ivar),first_index(ivar),:);
                end
                ar_data=[ar_data,squeeze(this_auto)];
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
            if ii==1 && nmod
                legend(modelsLegs)
            end
        end
        orient(fig,'tall')
        myfigure=struct('name',fig,'title',fig_title,'scale',0.85);
        figure(xrep,myfigure);
        xrep.pagebreak();
    end
end
%% Smoothed shocks

if add_smoothed_shocks && ~isempty(obj.filtering)
    add_str='';
    for ifig=1:n_myshocks_figs
        if n_myshocks_figs>1
            add_str=['(',int2str(ifig),')'];
        end
        fig_title=['Smoothed shocks',add_str];
        fig=figure('name',fig_title);
        [~,r,c]=number_of_rows_and_columns_in_figure(ifig,nshocks,graph_nrows,graph_ncols);
        for ivar=(ifig-1)*graph_nstar+1:min(nshocks,ifig*graph_nstar)
            ii=ivar-(ifig-1)*graph_nstar;
            vshock=shock_list{ii};
            myshock_data=obj(1).filtering.Expected_smoothed_shocks.(vshock);
            for imod=2:numel(obj)
                myshock_data=[myshock_data,obj(imod).filtering.Expected_smoothed_shocks.(vshock)];
            end
            subplot(r,c,ii)
            plot(myshock_data,'linewidth',2)
            axis tight
            title(shock_list_tex_names{ivar},...
                'fontsize',12)
            if ii==1
                leg=legend({obj.filename});
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
if add_shocks_correlations && ~isempty(obj.filtering)
    for imod=1:numel(obj)
        shock_corr=corrcoef(rise_time_series.collect(obj(imod).filtering.Expected_smoothed_shocks));
        mytable=[
            [' ',shock_list]
            shock_list',num2cell(shock_corr)
            ];
        mytable=struct('title',[obj(imod).legend,' model: Shock correlation structure '],...
            'table',{mytable});
        table(xrep,mytable);
        xrep.pagebreak();
    end
end
%% empirical distribution of shocks

if add_empirical_distribution_of_shocks && ~isempty(obj.filtering)
    add_str='';
    for ifig=1:n_myshocks_figs
        if n_myshocks_figs>1
            add_str=['(',int2str(ifig),')'];
        end
        fig_title=['Empirical PDF of shocks',add_str];
        fig=figure('name',fig_title);
        [~,r,c]=number_of_rows_and_columns_in_figure(ifig,nshocks,graph_nrows,graph_ncols);
        for ivar=(ifig-1)*graph_nstar+1:min(nshocks,ifig*graph_nstar)
            ii=ivar-(ifig-1)*graph_nstar;
            vshock=shock_list{ii};
            myshock_data=cell(1,numel(obj));
            lb=inf;
            ub=-inf;
            for imod=1:numel(obj)
                myshock_data{imod}=double(obj(imod).filtering.Expected_smoothed_shocks.(vshock));
                lb=min(lb,min(myshock_data{imod}));
                ub=max(ub,max(myshock_data{imod}));
            end
            myff=[];
            for imod=1:numel(obj)
                [ff_,xx_]=distributions.kernel_density(myshock_data{imod},lb,ub,'normal',250);
                myff=[myff,ff_(:)];
            end
            subplot(r,c,ii)
            plot(xx_(:),myff,'linewidth',2)
            axis tight
            title(shock_list_tex_names{ivar},'fontsize',12)
            if ii==min(nshocks,ifig*graph_nstar)
                leg=legend(modelsLegs);%,'location','SO','orientation','horizontal'
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

