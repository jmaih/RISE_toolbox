function [data,p]=create_data(p)
if nargin==0
    p=struct();
end

gpr=load('gpr.dat');

p.zss = exp(mean(gpr(:,1)));
p.paiss = exp(mean(gpr(:,2)));
p.beta = (p.zss*p.paiss)/(exp(mean(gpr(:,3))));

start_date='1983Q1';
gt = gpr(:,1) - log(p.zss);
pit = gpr(:,2) - log(p.paiss);
rt = gpr(:,3) - log(p.zss) + log(p.beta) - log(p.paiss);

data=ts(start_date,[gt,pit,rt],{'GHAT','PAIHAT','RHAT'});
data=pages2struct(data);


