%% housekeeping

% look up help for rfvar/structural_form

% Warning: With too many restrictions, there might not be a suitable
% rotation leading to a solution

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
%% Set the RFVAR
clc
tmpl=rfvar.template();
tmpl.constant=true;
tmpl.nlags=3;
tmpl.endogenous={'robs','"interest rates"','dy','"GDP growth"',...
    'labobs','"hours worked"','pinfobs','"inflation"',...
    'dw','"wages growth"','dc','dinve'};
%%
rv=rfvar(tmpl,'data',data);

%% Estimate the reduced-form VAR

rv1=estimate(rv);

%% add some options

structural_shocks={'mp','"Monetary Policy"','ad','"Aggregate demand"',...
    'as','"Aggregate supply"','ls','"Labor supply"'};
restrict_irf_zero={
    %     'dy{inf}@mp',0
    %     'dy{inf}@ad',0
    %     'dy{inf}@as',0
    %     'dy{inf}@ls',0
    };
restrict_irf_sign={
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

rv1=set(rv1,'structural_shocks',structural_shocks,...
    'restrict_irf_zero',restrict_irf_zero,...
    'restrict_irf_sign',restrict_irf_sign);

% note we could have added restrictions on the lag structure. In that case
% the syntaxt would be 'restrict_lags',{...}

%% Check identification
rv1=check_identification(rv1);

%% set ONE structural form

rv11=structural_form(rv1);

%% Plot the IRFs for one particular rotation

myirfs=irf(rv11);

sstate=get(rv11,'sstate');

shock_names=rv11.exogenous.name;
shock_texnames=rv11.exogenous.tex_name;
var_names=rv11.endogenous.name;
var_texname=rv11.endogenous.tex_name;
for ishock=1:numel(shock_names)
    shock=shock_names{ishock};
    figure('name',['Orthogonalized responses to a ',...
        shock_texnames{ishock},' shock']);
    for ivar=1:numel(var_names)
        vname=var_names{ivar};
        m=sstate.(vname);
        if m
            s=m;
        else
            s=1;
        end
        subplot(3,3,ivar)
        plot((myirfs.(shock).(vname)-m)/s)
        title(var_texname{ivar})
    end
end


