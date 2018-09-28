function xrep1=Tao_report(m,rep_data,my_varlist,unobsList)
%   Tao_report example of report
%  
%   Syntax
%   -------
%   ::
%       Tao_report(m,rep_data)
%  
%       Tao_report(m,rep_data,my_varlist)
%  
%       Tao_report(m,rep_data,my_varlist,unobsList)
%  
%   Inputs
%   -------
%  
%  - m :[dsge object|rise object] could be one or several objects
%  
%  - rep_data : [struct] with fields % 'title','address','author','email','abstract'
%       - name (optional) [char]: name under which the report will be saved
%       - title (optional) [char|cellstr]: title of the report
%       - address (optional) [char|cellstr]: address of the author(s)
%       - author (optional) [char|cellstr]: name(s) of the author(s)
%       - email (optional) [char|cellstr]: email(s) of the author(s)
%       - abstract (optional) [char|cellstr]: abstract of the report
%  
%  - my_varlist : [char|cellstr] : list of the variables of interest
%  
%  - unobsList : [char|cellstr] : list of the unobservable variables of interest
%  
%   Outputs
%   --------
%
%   - xrep1 : handle to the report batch. 
%           N.B: if requested, no report will be published
%  
%   More About
%   ------------
%  
%   Examples
%   ---------
%  
%   See also: 

obsList=m(1).observables.name;
if nargin<4
    unobsList=setdiff(m(1).endogenous.name,obsList);
    check_default_variables(m(1),unobsList,'unobservable')
    if nargin<3
        my_varlist=m(1).endogenous.name;
        check_default_variables(m(1),my_varlist,'endogenous')
        if nargin<2
            error('You need to provide at least two arguments')
        end
    end
end
close all

filtration=filter(m);

%% Make some choices

do_model=true;
do_model_equations=true;
do_steady_state_results=true;
do_estimation_results=true;
do_data_plot=true;
do_data_against_transition_probabilities=true;
do_generalized_irfs=true;
do_simple_irfs=true;
do_historical_decomposition=true;
do_counterfactual=~true;
do_variance_decomposition=true;
do_forecast=true;
do_shock_correlations=true;
do_empirical_distribution_of_shocks=true;
do_smoothed_shocks=true;
do_vector_autocorrelations=true;
do_smoothed_probabilities=true;
do_plot_unobservables=true;
graph_nrows=4;
graph_ncols=4;

irf_shocks=m(1).exogenous.name;

number_of_models=numel(m);
model_tags=cell(1,number_of_models);
for imod=1:number_of_models
    model_tags{imod}=m(imod).legend;
    if isempty(model_tags{imod})
        model_tags{imod}=m(imod).filename;
        % model_tags{imod}=sprintf('model_%0.0f',imod);
    end
end
earliest_date=m(1).options.estim_start_date;
latest_date=m(1).options.estim_end_date;

%% start the report with a title page
xrep=newreport('date',datestr(now));
myfields={'name','title','address','author','email','abstract'};
for ifield=1:numel(myfields)
    ff=myfields{ifield};
    if isfield(rep_data,ff)
        xrep.(ff)=rep_data.(ff);
    end
end

xrep.pagebreak()

%% locate the unobservables

unobs_locs=locate_variables(unobsList,m(1).endogenous.name);
unobsList_tex_names=get_tex_names(unobsList,m(1).endogenous.tex_name(unobs_locs));

%% locate the variables of interest
locs=locate_variables(my_varlist,m(1).endogenous.name);
var_list_tex_names=get_tex_names(my_varlist,m(1).endogenous.tex_name(locs));
graph_nstar=graph_nrows*graph_ncols;
nvar=numel(my_varlist);

%% get the list of the observed variables
obsListTexnames=get_tex_names(obsList,m(1).observables.tex_name);

nvarobs=numel(obsList);
n_myobs_figs=ceil(nvarobs/graph_nstar);

%% exogenous
shock_list=m(1).exogenous.name;
shock_list_tex_names=get_tex_names(shock_list,m(1).exogenous.tex_name);
%% add the model equations
if do_model
    xrep.section('title','Model code')
    report(m(end),xrep,'code')
    xrep.pagebreak()
end
%% do various descriptions
xrep.section('title','Description of variables')
report(m(1),xrep,'endogenous')
report(m(1),xrep,'exogenous')
report(m(1),xrep,'observables')
xrep.pagebreak()
%% add the model equations
if do_model_equations
    xrep.section('title','Model equations')
    report(m(1),xrep,'equations')
end
%% add the estimation results
if do_estimation_results && ~isempty(m(1).estimation)
    xrep.section('title','Estimation results')
    report(m,xrep,'estimation')
    xrep.pagebreak()
    report(m,xrep,'estimation_statistics')
end
%% add the steady state results
if do_steady_state_results
    steady_table=cell(nvar+2,1);
    tmp=struct2cell(var_list_tex_names);
    steady_table(3:end,1)=tmp(:);
    big_table=cell(1,numel(m)+1);
    big_table{1}=steady_table;
    for imod=1:numel(m)
        sstate=get(m(imod),'sstate');
        regs_i=m(imod).markov_chains.regime_tex_names;
        nregs=numel(regs_i);
        tble=cell(nvar+2,nregs);
        tble{1,1}=model_tags{imod};
        tble(2,1:nregs)=regs_i;
        for ivar=1:nvar
            tble(ivar+2,1:nregs)=num2cell(sstate.(my_varlist{ivar}));
        end
        big_table{1+imod}=tble;
    end
    xrep.table('title','Steady state values','log',[big_table{:}]);
    xrep.pagebreak()
end

%% plot the data

if do_data_plot
    db=pages2struct(m(1).options.data);
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
            plot([earliest_date,':',latest_date],db.(obsList{ivar}),'linewidth',2);
            title(obsListTexnames.(obsList{ivar}));
        end
        xrotate(45)
        xrep.figure('name',fig,'title',fig_title)
        xrep.pagebreak()
    end
end

%% plot the smoothed probabilities for all models
if do_smoothed_probabilities
    for imod=1:number_of_models
        [state_names,state_tex_names,nstates]=get_state_names(m(imod));
        fig_title=sprintf('Smoothed probabilities for %s model',model_tags{imod});
        fig=figure('name',fig_title);
        [rr,cc]=utils.plot.one_figure_rows_columns(nstates);
        for ivar=1:nstates
            subplot(rr,cc,ivar)
            vname=state_names{ivar};
            plot(filtration{imod}.smoothed_state_probabilities.(vname),'linewidth',2)
            ylim([-0.01,1.01])
            title(sprintf('%s(%s)',state_tex_names{ivar},vname))
        end
        xrotate(45)
        xrep.figure('title',fig_title,'name',fig)
        xrep.pagebreak()
    end
end

%% plot the observables against the transition probabilities
if do_data_against_transition_probabilities
    for imod=1:number_of_models
        smoothed=filtration{imod}.Expected_smoothed_variables;
        [state_names,state_tex_names,nstates]=get_state_names(m(imod));
        smooth_probs=filtration{imod}.smoothed_state_probabilities;
        for icase=1:nstates
            titel=sprintf('%s:: Observed vs %s(%s)',...
                model_tags{imod},state_tex_names{icase},state_names{icase});
            fig_handles=utils.plot.multiple(...
                @(xname)one_observable_against_transition(xname,obsListTexnames,smoothed,smooth_probs,state_names{icase}),...
                obsList,titel,graph_nrows,graph_ncols,...
                'FontSize',11,'FontWeight','normal');
            record_figures(xrep,fig_handles,titel,45);
        end
    end
end
%% plot unobserved variables
if do_plot_unobservables
    titel='Unobserved variables';
    fig_handles=utils.plot.multiple(...
        @(xname)one_unobservable(xname,unobsList_tex_names,m,model_tags),...
        unobsList,titel,graph_nrows,graph_ncols,...
        'FontSize',11,'FontWeight','normal');
    record_figures(xrep,fig_handles,titel,45);
end
%% compute impulse responses
%  list of shocks and preferred variables
if do_generalized_irfs
    myirfs=irf(m,'irf_periods',40,'irf_type','girf',...
        'irf_regime_specific',false);
    for ishock=1:numel(shock_list)
        shock=shock_list{ishock};
        if ~ismember(shock,irf_shocks)
            continue
        end
        titel=sprintf('(Generalized) IRFs to a %s shock',shock_list_tex_names.(shock));
        fig_handles=utils.plot.multiple(...
            @(xname)one_irf(xname,var_list_tex_names,myirfs.(shock),model_tags),...
            my_varlist,titel,graph_nrows,graph_ncols,...
            'FontSize',11,'FontWeight','normal');
        record_figures(xrep,fig_handles,titel);
    end
end

%% simple regime-specific irfs
if do_simple_irfs
    reg_names=m(1).markov_chains.regime_names;
    do_regime_specific=true;
    for imod=2:numel(m)
        do_regime_specific=do_regime_specific &&...
            isequal(reg_names,m(imod).markov_chains.regime_names);
    end
    
    do_regime_specific=do_regime_specific && numel(m)>1 && numel(reg_names)>1;
    
    if do_regime_specific
        % describe the regimes
        %---------------------
        reg_table=m(1).markov_chains.regimes;
        for ireg=1:m(1).markov_chains.regimes_number
            for ichain=1:m(1).markov_chains.chains_number
                chain_name=sprintf('%s_%0.0f',reg_table{1,ichain+1},reg_table{ireg+1,ichain+1});
                loc=strcmp(chain_name,m(1).markov_chains.state_names);
                chain_name_tex=m(1).markov_chains.state_tex_names{loc};
                reg_table{ireg+1,ichain+1}=chain_name_tex;
            end
        end
        xrep.table('title','Description of regimes',...
            'log',reg_table);
        xrep.pagebreak()
        
        % do the irfs
        %------------
        myirfs=irf(m,'irf_periods',40);
        for ireg=1:numel(reg_names)
            this_reg_name=reg_names{ireg};
            for ishock=1:numel(shock_list)
                shock=shock_list{ishock};
                if ~ismember(shock,irf_shocks)
                    continue
                end
                titel=sprintf('%s:: IRFs to a %s shock',this_reg_name,shock_list_tex_names.(shock_list{ishock}));
                fig_handles=utils.plot.multiple(...
                    @(xname)one_irf(xname,var_list_tex_names,myirfs.(shock).(this_reg_name),model_tags),...
                    my_varlist,titel,graph_nrows,graph_ncols,...
                    'FontSize',11,'FontWeight','normal');
                record_figures(xrep,fig_handles,titel);
            end
        end
    end
end
%% historical decomposition of shocks
if do_historical_decomposition
    histdec=historical_decomposition(m);
    [~,r,c]=number_of_rows_and_columns_in_figure(1,number_of_models,graph_nrows,graph_ncols);
    if ~iscell(histdec)
        histdec={histdec};
    end
    for ivar=1:nvar
        vname=my_varlist{ivar};
        if ivar==1
            contrib_names=histdec{1}.(vname).varnames;
            locs=locate_variables(contrib_names,shock_list,true);
            for ii=1:numel(locs)
                if isnan(locs(ii))
                    continue
                end
                contrib_names{ii}=shock_list_tex_names.(shock_list{locs(ii)});
            end
        end
        fig_title=['historical decomposition of ',var_list_tex_names.(my_varlist{ivar})];
        fig=figure('name',fig_title);
        for imod=1:number_of_models
            subplot(r,c,imod)
            plot_decomp(histdec{imod}.(vname))
            title(model_tags{imod},'interpreter','none')
        end
        xrotate(45)
        legend(contrib_names)
        xrep.figure('name',fig,'title',fig_title)
        xrep.pagebreak()
    end
end

%% counterfactual: what if the economy had been in the low/high volatility regime all the time
if do_counterfactual
    for imod=1:number_of_models
        % actual,counterf need to be defined
        titel=[model_tags{imod},'::Counterfactual: low volatility and random-walk regime'];
        fig_handles=utils.plot.multiple(...
            @(xname)one_counterfactual(xname,var_list_tex_names,actual,counterf),...
            my_varlist,titel,graph_nrows,graph_ncols,...
            'FontSize',11,'FontWeight','normal');
        record_figures(xrep,fig_handles,titel);
    end
end
%% variance decomposition
if do_variance_decomposition
    vardec=variance_decomposition(m);
    if ~iscell(vardec)
        vardec={vardec};
    end
    [r,c]=utils.plot.one_figure_rows_columns(number_of_models);
    % plot the decomposition
    for ivar=1:nvarobs
        vname=obsList{ivar};
        fig_title=['Variance decomposition of ',obsListTexnames.(obsList{ivar})];
        fig=figure('name',fig_title);
        for imod=1:number_of_models
            subplot(r,c,imod)
            tmp=100*vardec{imod}.conditional.(vname);
            tmp.varnames=vardec{imod}.conditional.(vname).varnames;
            % plt=plot_decomp(tmp(1:24)); did not work and did not return
            % an error. This begs some investigating
            plt=plot_decomp(tmp('1:24'));
            title([obsListTexnames.(obsList{ivar}),'(',model_tags{imod},')'],'interp','none')
            if imod==1
                contrib_names=vardec{imod}.conditional.(vname).varnames;
                locs=locate_variables(contrib_names,shock_list,true);
                for ii=1:numel(locs)
                    if isnan(locs(ii))
                        continue
                    end
                    contrib_names{ii}=shock_list_tex_names.(shock_list{locs(ii)});
                end
                legend(contrib_names,'location','NorthEast');
            end
        end
        xrep.figure('name',fig,'title',fig_title)
        xrep.pagebreak()
    end
end
%% real-time forecasting performance of the estimated model
if do_forecast
    [ts_fkst,ts_rmse]=forecast_real_time(m);
    if ~iscell(ts_fkst),ts_fkst={ts_fkst}; end
    pp=[];
    sstates=cell(1,number_of_models);
    ergodic_probs=cell(1,number_of_models);
    [r,c]=utils.plot.one_figure_rows_columns(number_of_models+1);
    for ivar=1:nvarobs
        fig_title=['real-time forecasts for ',obsListTexnames.(obsList{ivar})];
        fig=figure('name',fig_title);
        for imod=1:number_of_models
            if ivar==1
                sstates{imod}=cell2mat(m(imod).solution.ss);
                Qi=m(imod).solution.transition_matrices.Q;
                nreg=size(Qi,1);
                ergodic_probs{imod}=[eye(nreg)-Qi';ones(1,nreg)]\[zeros(nreg,1);1];
            end
            subplot(r,c,imod)
            [h]=plot_real_time(ts_fkst{imod}.(obsList{ivar}),pp);
            % add the steady state
            vloc=locate_variables(obsList{ivar},m(imod).endogenous.name);
            hold on
            xlim=get(h,'xlim');
            ergodic_mean=sstates{imod}(vloc,:);
            if numel(ergodic_mean)>1
                ergodic_mean=ergodic_mean*ergodic_probs{imod};
            end
            plot(xlim,ergodic_mean*ones(1,2),'color',[1,0,0])
            hold off
            title([obsListTexnames.(obsList{ivar}),'(',model_tags{imod},')'],'interpreter','none')
            grid on
        end
        xrotate(45)
        % plot the rmses
        subplot(r,c,number_of_models+1)
        plot(ts_rmse.(obsList{ivar}))
        title('RMSE')
        leg=legend(model_tags);
        set(leg,'interpreter','none')
        xrep.figure('name',fig,'title',fig_title)
        xrep.pagebreak()
    end
end
%% vector autocorrelations
if do_vector_autocorrelations
    [Auto]=theoretical_autocorrelations(m,'autocorr_ar',40);
    if ~iscell(Auto)
        Auto={Auto};
    end
    locs=locate_variables(my_varlist,m(1).endogenous.name);
    indices=get_tex_names(my_varlist,num2cell(locs));
    ar=size(Auto{1},3);
    xx=(1:ar)';
    mynames=cell(nvar^2,1);
    iter=0;
    for ivar=1:nvar
        for jvar=1:nvar
            iter=iter+1;
            mynames{iter}=[my_varlist{ivar},'@',my_varlist{jvar}];
        end
    end
    titel='Vector autocorrelations';
    fig_handles=utils.plot.multiple(...
        @(xname)one_vector_autocorrelation(xname,var_list_tex_names,indices,Auto,xx,model_tags),...
        mynames,titel,graph_nrows,graph_ncols,...
        'FontSize',11,'FontWeight','normal');
    record_figures(xrep,fig_handles,titel);
end

%% smoothed shocks

if do_smoothed_shocks
    titel='Smoothed shocks';
    fig_handles=utils.plot.multiple(...
        @(xname)one_smoothed_shock(xname,shock_list_tex_names,m,model_tags),...
        shock_list,titel,graph_nrows,graph_ncols,...
        'FontSize',11,'FontWeight','normal');
    record_figures(xrep,fig_handles,titel,45);
end
%% correlation of shocks

if do_shock_correlations
    report(m(1),xrep,'exogenous')
    xrep.pagebreak()
    for imod=1:number_of_models
        tmp=ts.collect(filtration{imod}.Expected_smoothed_shocks);
        shock_corr=corrcoef(tmp);
        shock_locs=locate_variables(shock_list,tmp.varnames);
        mytable=[
            [' ',shock_list(:).']
            shock_list(:),num2cell(shock_corr(shock_locs,shock_locs))
            ];
        xrep.table('title',['Shock correlation structure in ',model_tags{imod}],...
            'log',mytable);
        xrep.pagebreak()
    end
end

%% empirical distribution of shocks

if do_empirical_distribution_of_shocks    
    titel='Empirical distribution of smoothed shocks';
    fig_handles=utils.plot.multiple(...
        @(xname)one_smoothed_shock_pdf(xname,shock_list_tex_names,m,model_tags),...
        shock_list,titel,graph_nrows,graph_ncols,...
        'FontSize',11,'FontWeight','normal');
    record_figures(xrep,fig_handles,titel);
end
%% publish the report
if nargout
    xrep1=xrep;
else
    publish(xrep);
    close all
end
end

function check_default_variables(m,unobsList,type)
    tmp=locate_variables(unobsList,m.endogenous.name,'silent');
    disp(unobsList)
    if any(isnan(tmp))
        error(['A list of ',type,' has not been specified and some of the variables(if not all) in my default list above are not among the endogenous variables'])
    else
        warning('You have not specified a list of unobservable variables. I''ll fall back on the default list displayed above')
        pause(3)
    end
end

function record_figures(xrep,fig_handles,titel,angle)
if nargin<4
    angle=0;
end
add_str='';
nfigs=numel(fig_handles);
for ifig=1:nfigs
    if nfigs>1
        add_str=['(',int2str(ifig),')'];
    end
    fig_title=sprintf('%s%s',titel,add_str);
    fig=fig_handles(ifig);
    figure(fig)
    if angle
        xrotate(angle)
    end
    xrep.figure('name',fig,'title',fig_title)
    xrep.pagebreak()
end
end

function txnames=get_tex_names(names,tex_names)
txnames=struct();
for iname=1:numel(names)
    txnames.(names{iname})=tex_names{iname};
end
end

function [tex_name,legend_]=one_smoothed_shock_pdf(vshock,shock_list_tex_names,m,model_tags)
number_of_models=numel(m);
myshock_data=cell(1,number_of_models);
lb=inf;
ub=-inf;
for imod=1:number_of_models
    myshock_data{imod}=double(filtration{imod}.Expected_smoothed_shocks.(vshock));
    lb=min(lb,min(myshock_data{imod}));
    ub=max(ub,max(myshock_data{imod}));
end
for imod=1:number_of_models
    [ff_,xx_]=distributions.kernel_density(myshock_data{imod},lb,ub,'normal',250);
    if imod==1
        ff_=ff_(:);
        myff=ff_(:,ones(1,number_of_models));
    else
        myff(:,imod)=ff_(:);
    end
end
plot(xx_(:),myff,'linewidth',2)
axis tight
tex_name=shock_list_tex_names.(vshock);
legend_=model_tags;
end

function [tex_name,legend_]=one_smoothed_shock(vshock,shock_list_tex_names,m,model_tags)
myshock_data=filtration{1}.Expected_smoothed_shocks.(vshock);
for imod=2:numel(m)
    myshock_data=[myshock_data,filtration{imod}.Expected_smoothed_shocks.(vshock)]; %#ok<AGROW>
end
plot(myshock_data,'linewidth',2)
axis tight
tex_name=shock_list_tex_names.(vshock);
legend_=model_tags;
end

function [tex_name,legend_]=one_vector_autocorrelation(vname,var_list_tex_names,indices,Auto,xx,model_tags)
arobase=find(vname=='@');
vname_1=vname(1:arobase-1);
vname_2=vname(arobase+1:end);
first_index=indices.(vname_1);
second_index=indices.(vname_2);
ar_data=[];
for imod=1:numel(Auto)
    ar_data=[ar_data,squeeze(Auto{imod}(first_index,second_index,:))]; %#ok<AGROW>
end
plot(xx,ar_data,'linewidth',2)
axis tight
ylabel(var_list_tex_names.(vname_1))
xlabel(sprintf('%s_{t-s}',var_list_tex_names.(vname_2)))
tex_name='';
legend_=model_tags;
end

function [tex_name,legend_]=one_counterfactual(vname,var_list_tex_names,actual,counterf)
plot([actual.(vname),counterf.(vname)],'linewidth',2)
tex_name=var_list_tex_names.(vname);
legend_={'actual','counterfactual'};
end

function [tex_name,legend_]=one_irf(vname,var_list_tex_names,irf_batch,model_tags)
plot(irf_batch.(vname),'linewidth',2)
tex_name=var_list_tex_names.(vname);
legend_=model_tags;
end

function [tex_name,legend_]=one_unobservable(vname,unobsList_tex_names,m,model_tags)
V=filtration{1}.Expected_smoothed_variables.(vname);
for imod=2:numel(m)
    V=[V,filtration{imod}.Expected_smoothed_variables.(vname)];
end
plot(V,'linewidth',2)
tex_name=unobsList_tex_names.(vname);
legend_=model_tags;
end

function [tex_name,legend_]=one_observable_against_transition(vname,all_tex_names,smoothed,smooth_probs,which_state)
dd=smoothed.(vname);
[ax,h1,h2]=plotyy(dd,smooth_probs.(which_state),'linewidth',2);
set(ax(1),'YLim',[min(dd) max(dd)])
set(ax(2),'YLim',[-0.01 1.01])
set(ax(1),'Box','off')
tex_name=all_tex_names.(vname);
legend_={'data','prob'};
end

function [state_names,state_tex_names,nstates]=get_state_names(model)

state_names=model.markov_chains.state_names;
state_tex_names=model.markov_chains.state_tex_names;
const_loc=strcmp(state_names,'const_1');
state_names(const_loc)=[];
state_tex_names(const_loc)=[];
nstates=numel(state_names);
end