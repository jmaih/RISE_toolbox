%% load data and plot the data
%-----------------------------
data  = load('data/US_EURO.txt');
date_start    = '1983Q1';
vnames={'YGR','INFL','INT'};
Countries={'US','EURO'};
rawdata=struct();
for ic=1:numel(Countries)
    for iname=1:numel(vnames)
        rawdata.(Countries{ic}).(vnames{iname})=ts(date_start,data(:,(ic-1)*3+iname));
    end
end
%% plot/visualize the data
close all
country='all';
figure('name','data for the An and Schorfheide model');
for iname=1:numel(vnames)
    name=vnames{iname};
    legend_=country;
    switch lower(country)
        case 'us'
            dat_=rawdata.(Countries{1}).(name);
        case 'euro'
            dat_=rawdata.(Countries{2}).(name);
        case 'all'
            dat_=[rawdata.(Countries{1}).(name),rawdata.(Countries{2}).(name)];
            legend_=Countries;
        otherwise
            error(['unknown country ',country])
    end
    subplot(3,1,iname)
    plot(dat_)
    title(name)
    if iname==1
        legend(upper(legend_))
    end
end
%% RISE the various versions of model
m1=rise('as2007','rise_flags',{'gap_rule',true});
m2=rise('as2007','rise_flags',{'gap_rule',false});
%% variables of interest
myvars={'C','D','GBAR','AC','H','M','PAI','R','W','Y','YSTAR'};% B G LAMBDA N Q10 RSTAR SC T XI
%% do some impulse responses

M=set([m1,m2],'legend',{'gap rule','growth rule'});

myirfs=irf(M);
%% Plot irfs
nvars=numel(myvars);

for ishock=1:m1.exogenous.number(1)
    shock=m1.exogenous.name{ishock};
    tex_shock=m1.exogenous.tex_name{ishock};
    figure('name',['impulse responses to a ',tex_shock,' shock']);
    for ivar=1:nvars
        subplot(4,3,ivar)
        vname=myvars{ivar};
        vpos=strcmp(vname,m1.endogenous.name);
        tex_vname=m1.endogenous.tex_name{vpos};
        plot(myirfs.(shock).(vname),'linewidth',2)
        title(tex_vname)
        if ivar==1
            legend({M.legend})
        end
    end    
end
%% simulate artificial data
close all
mysimul=simulate(M);
figure('name','Simulated data');
for ivar=1:nvars
    subplot(4,3,ivar)
    vname=myvars{ivar};
    vpos=strcmp(vname,m1.endogenous.name);
    tex_vname=m1.endogenous.tex_name{vpos};
    plot(mysimul.(vname),'linewidth',2)
    title(tex_vname)
    if ivar==1
        legend({M.legend})
    end
end
%% variance decomposition
myvdec=variance_decomposition(M);

for imod=1:numel(M)
    figure('name',['Variance Decomposition in ',M(imod).legend]);
    for ivar=1:nvars
        subplot(4,3,ivar)
        vname=myvars{ivar};
        vpos=strcmp(vname,m1.endogenous.name);
        tex_vname=m1.endogenous.tex_name{vpos};
        plot_decomp('1:30',myvdec{imod}.conditional.(vname))
        title(tex_vname)
    end
end
%% estimate a model
%-------------------
clc
m1est=estimate(m1,'data',rawdata.US);

%% historical decomposition
mydec=historical_decomposition(m1est)
