%% housekeeping
close all
clear
clc
%% Information for the title page

% N.B: Be careful with the "_" it has to be escaped with a "\" in order for
% latex to process it as we wish.

rep_data=struct();
myfields={'name','title','address','author','email','abstract'};

rep_data.name='newkeynesswitch'; % Name under which the report will be saved
rep_data.title='Switches in the US Macroeconomic Data: Policy or Volatility?';
rep_data.address='Your institution'; 
rep_data.author='first\_name last\_name'; % your name 
rep_data.email='first\_name@last\_name.com'; % your email
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
rep_data.abstract=myabstract;
%% load the estimated models and relabel them if necessary
nk_models={'volatilityOnly','policyOnly','volatilityPolicySame','volatilityPolicyIndependent'};
nk_models_newlabels={'volOnly','polOnly','volPolSame','volPolInd'};

model_objects=rise.empty(0);
for imod=1:numel(nk_models)
    tmp=load([nk_models{imod},filesep,'estimation',filesep,'estimated_model']);
    obj=tmp.obj;
    model_objects(imod,1)=obj;
    model_objects(imod,1).legend=nk_models_newlabels{imod};
end
Tao_report(model_objects,rep_data)

