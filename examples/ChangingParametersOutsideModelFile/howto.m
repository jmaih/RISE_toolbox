clc
clear classes
disp(upper('this file may need updating... contact junior.maih@gmail.com if it does not work'))

test=rise('rwz12',...
    'steady_state_file','rwz_steady_state',...
    'irf_anticipate',1,'irf_type','girf');
%     ,...
%%
clear m
N=5;
rhor_vals=linspace(0,.99,N);
m(1:N)=test;
for ii=1:N
    pp=struct('name','rhor','value',rhor_vals(ii),'regime',2);
    m(ii)=m(ii).set_parameters(pp);
end
%%

m=m.solve;
%%
m=m.irf;
%%
m(1)=m(1).irf;

%%
m.print_solution
