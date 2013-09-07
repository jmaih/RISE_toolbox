%% housekeeping
clear all
close all
clc
%% Instructions
% - Please run this file block by block and make sure you read the
% comments in each block to understand what it does. If there is anything
% you do not understand, ask questions to the instructor or to your
% neighbor
% - to run a particular block, click on the block and then on your keyboard
% press CTRL+Enter
%% familiarize yourself with RISE.
% - issue the command: methods(rise)
% - create an empty rise object: junk=rise.empty(0)
% - read the different field names of junk and check that you understand
% the meaning of most of them
% - issue the command methods(junk)
% - call the following methods on the empty object (junk): irf,
% variance_decomposition, simulate, estimate, evaluate, solve, forecast,
% etc. do this by issuing e.g. the commands: irf(junk), solve(junk), etc.
% try to understand intuitively what rise does when you pass an empty
% object to a particular method.
%% Now, turn on the game
% open all the files with extension dyn in this folder and try to
% understand what they mean and how they are related to one another

%% Read the model(s)
% there are many parameterizations of the model and given the way the files
% are written we are going to take advantage of RISE's powerful
% macro-language, which allows conditional parsing of files.

% first we create some labels: the first column refers to indices that will
% be used when reading the models. the second one contains the actual
% labels
labels={0,'original'
    1,'policyChangeOnly'
    2,'volatilityChangeOnly'
    3,'privateChangeOnly'};

% we initialize an empty vector of models that will hold all models
m=rise.empty(0);

% we loop through the different models using the information in the labels
for imod=1:size(labels,1)
    m(imod,1)=rise('fwz',... % name of the file to read
        'rise_flags',struct('indx_model',labels{imod,1})... % conditional on this parameterization
        );
    % indx_model varies from 0 to 3 in the model file
end

%% We solve all the models simultaneously
m=solve(m);

%% check the stability of the system: remember there are 4 models

% refer to equation(22) in Dr Zha's slides

m.is_stable_system
 
% all the models in the m vector are stable. This is denoted by 1 (0 would
% mean unstable)
%% Print the solution of one model, say the first model in vector m

print_solution(m(1))

% to print the solution for one particular model, say the 3rd model in the
% m vector, write the following: print_solution(m(3))
%% Print the solution for all the models

print_solution(m)

%% Print the solution just for a subset of variables

print_solution(m,{'PAI','X','R'})

%% differences in solution approaches

for imod=1:numel(m)
    disp(['model ',labels{imod,2}])
    disp(max(max(max(abs(m(imod).T-mfwz(imod).T)))))
end

%% Compute impulse responses
myirfs=irf(m,... % all the models
    'irf_periods',15 ... % desired length of the irfs, default is 40
);

%% compare the irfs of the various models
shock_names=fieldnames(myirfs);
varlist={'PAI','R','X'};
for ireg=1:2
    regime=['regime_',sprintf('%0.f',ireg)];
    for ishock=1:numel(shock_names)
        shock=shock_names{ishock};
        figure('name',['impulse responses to a ',shock,' shock in ',regime])
        for ivar=1:numel(varlist)
            subplot(3,1,ivar)
            plot(myirfs.(shock).(regime).(varlist{ivar}),'linewidth',2)
            title(varlist{ivar},'interp','none')
            if ivar==1
                legend(labels(:,2))
            end
        end
        [~,tmp]=sup_label([shock,' shock, ',regime],'t');
        set(tmp,'fontsize',15)
    end
end
%% simulation of data
m=simulate(m);

%% recovering and plotting the simulations
% pick a model
choice=1;

% load all the simulations as a matrix
SIMULATIONS=vertcat(m(choice).varendo.value);

% recover the list of the endogenous variables from one model
varendo_list={m(1).varendo.name};

% choose a list of variables of interest. Let's just take the same list we
% had earlier.
newlist=varlist;

% plot the simulations for all those variables
figure('name',['simulated data for model ',sprintf('%0.f',choice)]);
for ivar=1:numel(newlist)
    subplot(3,1,ivar)
    % pick a variable
    v=newlist{ivar};
    % locate its position
    vloc=locate_variables(v,varendo_list);
    % plot it
    plot(SIMULATIONS(vloc,:),'linewidth',2)
    title(v,'interp','none')
end
[junk,tmp]=sup_label(['Simulated data in model ',sprintf('%0.f',choice)],'t');
set(tmp,'fontsize',15)
