%% housekeeping
clear all
close all
clc

%% initialize the report
instructions={};

%% add a title page
mytitle=struct('title','Asymmetric Expectation effects of regime shifts in monetary policy',...
    'address','San Francisco Fed & Atlanta Fed & Atlanta Fed ',...
    'date',datestr(now),...
    'author','Zheng Liu & Dan F. Waggoner & Tao Zha',...
    'email','zliu@imf.org & dwaggoner@imf.org & tzha@imf.org');
instructions=[instructions,{'title';mytitle}];

%% add some key points
myparagraph=struct('text',char('We have not succeeded in answering all of our problems',...
    'The answers we have found only serve to raise a whole set of new issues',...
    'In some ways we are as confused as ever. But we believe we are confused',...
    'on a higher level and about more important things'),...
    'title','Introduction',...
    'itemize',true);
instructions=[instructions,{'paragraph';myparagraph}];

%% add the model equations
instructions=[instructions,{'equations';[]}];
%% read the model file

lwz=rise('lwz09_2','irf_periods',20);

%% the naive economy: assume the Hawkish regime lasts forever
naive_calibration={'name','value','regime'
    'q_tp_2_1',0,1
    'q_tp_2_1',0,2
    };
lwz_naive=lwz.set_parameters(naive_calibration);
%% solve, evaluate, print solution, etc, not required

%% compute impulse responses
myirfs=irf(lwz);

%% compute irfs for the naive economy
my_naive_irfs=irf(lwz_naive);
%% plot impulse responses: naive vs normal

endo_names={lwz.varendo.name};
shock_names={lwz.varexo.name};
mylist={'PAI','Y','R'};
for ishock=1:numel(shock_names)
    shock=shock_names{ishock};
    fig_title=['IRFs to a ',lwz.varexo(ishock).tex_name,' shock (',shock,')'];
    theFig=figure('name',fig_title);
    for ivar=1:numel(mylist)
        vname=mylist{ivar};
        vloc=find(strcmp(vname,endo_names));
        subplot(3,1,ivar)
        expect=myirfs.(shock).(vname);
        naive=my_naive_irfs.(shock).(vname);
        tmp=plot([expect,naive],'linewidth',2);
        set(tmp(1:2),'linestyle','--')
        title(lwz.varendo(vloc).tex_name)
        if ivar==1
            legend({'dovish','hawkish','naive-dovish','naive-hawkish'})
        end
    end
    myfigure=struct('name',theFig,...
        'title',fig_title,'angle',90);
    instructions=[instructions,{'figure';myfigure}];
end

%% regime-dependent structural parameters ( eta and iota)

regDep_calibration={'name','value','regime'
    'eta',0.75,2
    'iota',0,2
    };

lwz=lwz.set_parameters(regDep_calibration);
lwz_naive=lwz_naive.set_parameters(regDep_calibration);

%% impulse responses
myirfs_2=irf(lwz);

my_naive_irfs_2=irf(lwz_naive);

%% plot IRFs
for ishock=1:numel(shock_names)
    shock=shock_names{ishock};
    fig_title=['IRFs (Regime Dep. structural params) to a ',lwz.varexo(ishock).tex_name,' shock (',shock,')'];
    theFig=figure('name',fig_title);
    for ivar=1:numel(mylist)
        vname=mylist{ivar};
        vloc=find(strcmp(vname,endo_names));
        subplot(3,1,ivar)
        expect=myirfs.(shock).(vname);
        naive=my_naive_irfs.(shock).(vname);
        tmp=plot([expect,naive],'linewidth',2);
        set(tmp(1:2),'linestyle','--')
        title(lwz.varendo(vloc).tex_name)
        if ivar==1
            legend({'dovish','hawkish','naive-dovish','naive-hawkish'})
        end
    end
    myfigure=struct('name',theFig,...
        'title',fig_title,'angle',90);
    instructions=[instructions,{'figure';myfigure}];
end

%%
Report_name='lwz_RED_2009';
retcode=rise_report(lwz,Report_name,instructions);
% % % retcode=rise_report(model,Report_name,instructions,graphicsPath);
