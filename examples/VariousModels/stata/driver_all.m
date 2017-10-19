%% housekeeping
clear
close all
clc

%%

noclearall=true;

txt=which('driver_all');

lastsep=find(txt==filesep,1,'last');

folder=txt(1:lastsep-1);

w=what(folder);

w=regexp(w.m,'driver\d{2}','match');

w=[w{:}];

bigm=rise.empty();

for ii=1:numel(w)
    
    run(w{ii})
    
    bigm(ii)=mest;
end

