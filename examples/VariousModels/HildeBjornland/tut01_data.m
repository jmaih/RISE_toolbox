%% housekeeping
clear
close all
clc % home

%% read the data
[d,t]=xlsread('bjornland_swe_Data.xlsx');

%% create variables and put them into the RISE time series format
vnames=t(6,2:end);

start='1982Q1';

db=ts(start,d,vnames);

db=pages2struct(db);

tex=struct();
tex.GDP='GDP';
tex.r='Domestic interest rate';
tex.CPISA='CPI seasonal adjusted';
tex.FED='Federal Funds rate';
tex.RERinv='Real effective exchange rate';
tex.Aldp='Annual inflation rate';
tex.du93Q1='Dummy';
tex.du95Q4='Dummy';
tex.du92Q3='Dummy';
tex.rtwi='Foreign interest rate';

%% create additional variables

db.LGDP=log(db.GDP);
tex.LGDP='Log GDP';

db.DGDP=db.LGDP-db.LGDP{-1};
tex.DGDP='GDP growth';

db.LRER=log(db.RERinv);
tex.LRER='Log RER';

db.DRER=db.LRER-db.LRER{-1};
tex.DRER='log change in RER';

%% plot the data (not the dummy variables)
mynames=fieldnames(db);

mynames=mynames-mynames(strncmp(mynames,'du9',3));

figure('name','database');
for ii=1:numel(mynames)
    v=mynames{ii};
    subplot(4,3,ii)
    plot(db.(v))
    title(tex.(v))
end
xrotate(45)

%% save data for later use
close all

save('tut01_data','db','tex')
