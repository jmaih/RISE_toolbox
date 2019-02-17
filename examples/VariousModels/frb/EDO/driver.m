%% housekeeping
clear
close all
clc

%% rise the model

m=rise('edo');

%% set parameters

[p,priors]=edo_params();

m=set(m,'parameters',p);

%% solve the model

m=solve(m,'steady_state_file','ssfile');

%% print_solution

print_solution(m)

%% do IRFs

myirfs=irf(m);

%% plot some irfs

% THIS NEEDS FIXING: QUICK IRFS DOES NO LONGER UNDERSTAND MULTIPLE
% VARIABLES IN ONE PLOT... TO BE FIXED
close all

myList={
    {'WC','WK'}%	"Detrended level of wages in the fast-growth sector"
    {'YC','YK'}%	"Detrended level of output in the fast-growth sector"
    'R'%	"Federal funds rate"
    {'HC','HSC'}
    'unemp'%	"Unemployment rate"
    'EIK'%	"Detrended level of business investment spending"
    'EC'%	"Detrended level of household spending on nondurables and nonhousing services"
    {'INFCNA','INFCOR'}%	"Core PCE inflation; PICXFE (/400) in FRB/US"
    {'GAP','PFGAP'}
    'RCH'%	"Rental rate for housing"
    'ECH'%	"Detrended level of household spending on residential construction"
    'QCH'%	"Market price of housing"
    };

shockList=[];

quick_irfs(set(m,'tex_name',{'EC','Detr. household spending on nondurables'
    'ECH','Detr. household spending on residential construction'}),...
    myirfs,myList,shockList,[4,4])