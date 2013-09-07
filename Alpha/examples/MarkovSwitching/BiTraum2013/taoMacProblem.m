% housekeeping
close all
clear all
clc

%% initialize the report with default options

xrep=rise_report('report_name','mac_problem');
xrep.pagebreak();

%% plot the log function

fig_title='log function';
fig=figure('name',fig_title);
x=1:100;
y=log(x);
plot(x,y,'linewidth',2)
myfigure=struct('name',fig,'title',fig_title,'scale',0.85);
figure(xrep,myfigure);
xrep.pagebreak();


%% now create the report
publish(xrep)
close all

