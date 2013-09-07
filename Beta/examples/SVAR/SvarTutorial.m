%% housekeeping
clear classes 
close all 
clc
%% load RISE
do_set_paths=true;
if do_set_paths
    addpath C:\Users\Junior\Documents\GitHub\RISE_toolbox
    rise_startup()
end
%% Load and transform data
testdata=peersman_data();

% Construct variables for use in VAR:
db=struct();
db.Dloil=100*log(testdata.oil_price/testdata.oil_price{-1});
db.Dlo=100*log(testdata.output_us/testdata.output_us{-1});
db.Dlp=100*log(testdata.cpi_us/testdata.cpi_us{-1});
db.r=testdata.intrate_us;
%% Set SVAR parameters
svar_order = 3; % Order of VAR lags

%% create a svar as a rise object
mysvar=struct('model','svar',...
    'var',{{'Dloil','$oil price inflation$','Dlo','$output growth$','Dlp','$CPI Inflation$','r','$Interest rate$'}});
% push the order of the VAR directly, and override the default of 4. You
% could always do it later as well
obj=rise(mysvar,'data',db,'svar_lags',svar_order);%
%% Cholesky identification
% by default, short-run Cholesky restrictions are applied in the order of
% endogenous variables as listed (alphabetically) in the model. Below, we
% impose a different order by specifying it explicitly.
obj=set_options(obj,'svar_order',{'Dloil','Dlp','Dlo','r'});
obj=set_properties(obj,'filename','short cholesky');
% here we could have explictly said we wanted a cholesky restriction
% obj=set_options(obj,'svar_restrictions','choleski'); but this is
% unnecessary since it is the default mode
%% Cholesky restrictions in the long run
% here we apply the same restrictions but in the long run
obj(2)=set_options(obj(1),'svar_restrictions','long_cholesky');
obj(2)=set_properties(obj(2),'filename','long cholesky');
%% replicating the short-run Cholesky identification through short-run restrictions
short_run_restrictions={ % eqtn,shock,restr_type,value
    'Dloil','Dlo','short',0
    'Dloil','Dlp','short',0
    'Dloil','r','short',0
    'Dlp','Dlo','short',0
    'Dlp','r','short',0
    'Dlo','r','short',0
    };
% push those restrictions into a 3rd object based on the first
% in this case, the option 'svar_order' from the previous model will not
% matter at all.
obj(3)=set_options(obj(1),'svar_restrictions',short_run_restrictions);
obj(3)=set_properties(obj(3),'filename','short run(cholesky)');
%% construct some exclusion restrictions, mixing short and long run
mixed_restrictions={ % eqtn,shock,restr_type,value
    'Dloil','Dlo','short',0
    'Dloil','Dlp','short',0
    'Dloil','r','short',0
    'Dlo','r','short',0
    'Dlo','r','long',0
    'Dlo','Dlo','long',0
    };
% push those restrictions into a 4th object based on the first
obj(4)=set_options(obj(1),'svar_restrictions',mixed_restrictions); % could also use short_run
obj(4)=set_properties(obj(4),'filename','mixed restrictions');
%% Solve all 4 models simultaneously
clc
objs=solve(obj,'debug',true,'optimset',struct('Display','off'),'svar_restarts',1);
%% compute irfs for all 3 models simultaneously
testirfs=irf(obj,'irf_periods',40);
%% plot the irfs
figure('name','impulse response functions')
endo_names=fieldnames(testirfs.(obj(1).varexo(1).name));
iter=0;
for iname=1:numel(endo_names)
    v=endo_names{iname};
    vloc=find(strcmp(v,{obj(1).varendo.name}));
    vtex=obj(1).varendo(vloc).tex_name;
    for ishock=1:obj(1).NumberOfExogenous
        shock_name=obj(1).varexo(ishock).name;
        shock_texname=obj(1).varexo(ishock).tex_name;
        iter=iter+1;
        subplot(4,4,iter)
        tmp=testirfs.(shock_name).(v);
        if ~strcmp(v,'r')
            tmpdata=cumsum(double(tmp));
            tmp=rise_time_series(tmp.TimeInfo,tmpdata);
        end
        plot(tmp,'linewidth',2)
        if ismember(iter,(1:4:16))
            ylabel(vtex)
        end
        if ismember(iter,(1:4))
            title(shock_texname,'interpreter','none')
        end
        if iname==1 && ishock==1
            leg=legend({obj.filename});
        end
    end
end

