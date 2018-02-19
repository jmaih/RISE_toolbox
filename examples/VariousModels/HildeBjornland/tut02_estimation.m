%% housekeeping
close all
clear
clc
%% load the data

load tut01_data

%% set up VAR model

warning('This system is nonstationary because of LGDP !!!')

endog={'rtwi','LGDP','Aldp','r','LRER'};

nlags=4;

const=true;

exog={'du93Q1','du95Q4','du92Q3'};

v=rfvar(endog,exog,nlags,const); % formerly redfvar

%%
clc
ve=estimate(v,db,{db.LGDP.start,db.LGDP.finish});

%% should the Feds fund rate react to domestic variables?
linres={};

for ilag=1:nlags
    
    for iv=2:numel(endog)
        
        y=endog{iv};
        
        linres=[linres;{sprintf('b%0.0f(1,%s)=0',ilag,y)}];
        
    end
    
end

%%

ve_lr=estimate(v,db,{db.LGDP.start,db.LGDP.finish},[],linres);

%% save estimated models for later use

models=struct('ve',ve,'ve_lr',ve_lr);

save('tut02_estimation','models')
									
