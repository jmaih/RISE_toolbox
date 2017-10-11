%% housekeeping
clear
close all
clc

%%

noclearall=true;

w=what(pwd);

w=regexp(w.m,'driver\d{2}','match');

w=[w{:}];

bigm=rise.empty();

for ii=1:numel(w)
    
    run(w{ii})
    
    bigm(ii)=mest;
end

