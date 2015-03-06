function [data,p]=create_data(p)
if nargin==0
    p=struct();
end

gygprrp=load('gygprrp.dat');

p.zss = exp(mean(gygprrp(:,1)));
p.beta = exp(log(p.zss)-mean(gygprrp(:,3)));

start_date='1959Q1';

% per capital GDP growth
gyt = gygprrp(:,1) - log(p.zss);
% Inflation growth
gpt = gygprrp(:,2);
% Real interest rate
rrpt = gygprrp(:,3) - log(p.zss) + log(p.beta);

data=ts(start_date,[gyt,gpt,rrpt],{'GY_HAT','GPAI_HAT','RRPAI_HAT'});
data=pages2struct(data);

end