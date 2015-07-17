%% housekeeping
clear all
close all
clc
%% loading of the RISE toolbox
rise_startup()

%% Bring in the data: set them as time series
rawdata=load('data_SwitchingFed11_demeaned');
fields=fieldnames(rawdata);
init_date='1959Q1';
data=[];
for ii=1:numel(fields)
    data=[data,rawdata.(fields{ii})];
end
data=ts('1959Q1',data,fields);
%% Set the Svensson-Rudebusch SVAR
r=svar.template();
r.nlags=2;
r.constant=false;
r.endogenous={'Pie','"inflation"','Y','"GDP growth"','R','"interest rates"'};
%% chain will control the 3rd equation only
r.markov_chains(1)=struct('name','inter',...
    'states_expected_duration',[3+1i,3+1i],... % this implies we have two states
    'controled_parameters',{{'a0(3)','a1(3)','a2(3)','sig(3)'}});
% add restrictions
%-----------------
restrict_lags={
    % first equation
    coef(1,'R',0),0
    coef(1,'Y',0),0
    coef(1,'R',1),0
    coef(1,'Y',2),0
    coef(1,'R',2),0
    % second equation
    coef(2,'Pie',1),0
    coef(2,'Y',1),0
    coef(2,'Pie',2),0
    coef(2,'Y',2),0
    coef(2,'R',2),0
    };
% third equation
for istate=1:2
    restrict_lags=[restrict_lags
        {
        coef(3,'R',0,'inter',istate),0
        coef(3,'Pie',0,'inter',istate),0
        coef(3,'R',1,'inter',istate)+coef(3,'Pie',1,'inter',istate),0
        coef(3,'R',2,'inter',istate),0
        coef(3,'Pie',2,'inter',istate),0
        }
        ];
end
%% create the svar object and assign the restrictions and the data
rv=svar(r,'estim_linear_restrictions',restrict_lags,'data',data);

%% estimate the structural VAR model
rv1=estimate(rv);

%% alternative specifications
%{

%% chain will control the first equation only
r.markov_chains(1)=struct('name','inter',...
    'states_expected_duration',[3+1i,3+1i],...
    'controled_parameters',{{'a0(1)','a1(1)','a2(1)','sig(1)'}});
% add restrictions
%-----------------
restrict_lags=cell(0,2);
for istate=1:2
    restrict_lags=[restrict_lags
        {
        % first equation
        coef(1,'R',0,'inter',istate),0
        coef(1,'Y',0,'inter',istate),0
        coef(1,'R',1,'inter',istate),0
        coef(1,'Y',2,'inter',istate),0
        coef(1,'R',2,'inter',istate),0
        }
        ];
end
restrict_lags=[restrict_lags
    {
    % second equation
    coef(2,'Pie',1),0
    coef(2,'Y',1),0
    coef(2,'Pie',2),0
    coef(2,'Y',2),0
    coef(2,'R',2),0
    % third equation
    coef(3,'R',0),0
    coef(3,'Pie',0),0
    coef(3,'R',1)+coef(3,'Pie',1),0
    coef(3,'R',2),0
    coef(3,'Pie',2),0
    }];
rv(2)=svar(r,'estim_linear_restrictions',restrict_lags);
%% chain will control the variances only
r.markov_chains(1)=struct('name','inter',...
    'states_expected_duration',[3+1i,3+1i],...
    'controled_parameters',{{'sig'}});
% add restrictions
%-----------------
restrict_lags={
    % first equation
    coef(1,'R',0),0
    coef(1,'Y',0),0
    coef(1,'R',1),0
    coef(1,'Y',2),0
    coef(1,'R',2),0
    % second equation
    coef(2,'Pie',1),0
    coef(2,'Y',1),0
    coef(2,'Pie',2),0
    coef(2,'Y',2),0
    coef(2,'R',2),0
    % third equation
    coef(3,'R',0),0
    coef(3,'Pie',0),0
    coef(3,'R',1)+coef(3,'Pie',1),0
    coef(3,'R',2),0
    coef(3,'Pie',2),0
    };
rv(3)=svar(r,'estim_linear_restrictions',restrict_lags);
%}
