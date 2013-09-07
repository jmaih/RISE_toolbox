%% housekeeping
close all
clear all
clc
%% add the necessary paths
rise_startup()

%% get the data
[datt,navn]=xlsread('DataNEMO');
startdate=[int2str(datt(1,1)),'Q',int2str(datt(1,2))];
rawdata=rise_time_series(startdate,datt(:,3:end),char(navn(3:end)));
obs_names=char('DPQ_P_NW','D_GDP_NW','RN3M_NW');
cond_obs_names=char('RN3M_NW');
data=window(rawdata,startdate,'',obs_names);

%% construction of artificial conditional information 
% using a 3-variable VAR in 4 lags
lags=4;
nsteps=12;
Y=transpose(cell2mat(data.data(2:end,2:end)));
[nvar,smpl0]=size(Y);
yt=Y(:,lags+1:end);
X=[];
for ll=1:lags
    X=[X;Y(:,(lags+1:end)-ll)]; %#ok<AGROW>
end
% add a constant
smpl=size(yt,2);
X=[X;ones(1,smpl)];
A=yt/X;
% Compute the steady state
C=A(:,end);
ys=eye(nvar);
for ll=1:lags
    Al=A(:,(ll-1)*nvar+1:ll*nvar);
    ys=ys-Al;
end
ys=ys\C;
% Companion form
B=[A(:,1:end-1)
    eye(nvar*(lags-1)),zeros(nvar*(lags-1),nvar)];
% compute forecasts
big_ys=repmat(ys,lags,1);
yf=zeros(nvar,nsteps,smpl0-lags+1);
for jj=1:smpl0-lags+1
    y0=Y(:,jj-1+(1:lags));
    for t=1:nsteps
        y0=fliplr(y0);
        tmp=B*(y0(:)-big_ys);
        yf(:,t,jj)=ys+tmp(1:nvar);
        y0=[y0(:,2:end),yf(:,t,jj)];
    end
end

% put that information into some conditional database
% at the end of lags, we make forecasts for period lags+1. Thus the start
% date is lags
tmp=rise_date(startdate);
raw_cond_data=rise_time_series(tmp.observation_2_date(lags),permute(yf,[3,1,2]),data.varnames);
cond_data=raw_cond_data.window('','',cond_obs_names);

%% constant-parameter model
profile on
mc=rise('Canonical_Const','data',data);
mc=mc.estimate('optimizer',1);
profile off
profile viewer
%%
check_optimum(mc)
%% constant-parameter model with dsge_var with change in the optimization specs
mc_var=rise('Canonical_dsge_var','data',data,...
    'optimset',optimset('MaxNodes',20,'TolFun',sqrt(eps),...
    'MaxTime',3*60,'MaxIter',2000));
mc_var=mc_var.estimate('optimizer','bee_gate');
check_optimum(mc_var)

%% constant parameter model with conditioning information
% add conditional data, change the horizon to 4. In this case all shocks
% are anticipated 4 periods in advance
mc2=mc.set_options('cond_data_ct',cond_data,'solve_expect_order',4,...
    'optimset',optimset('MaxNodes',20,'TolFun',sqrt(eps),...
    'MaxTime',3*60,'MaxIter',2000));
mc2=mc2.estimate('optimizer','bee_gate');
check_optimum(mc2)

%% constant parameter model with conditioning information on shorter 
% horizon for the conditioning data
% add conditional data, change the horizon to 4. In this case all shocks
% are anticipated 4 periods in advance
mc3=mc.set_options('cond_data_ct',cond_data.window('','','',1:4),'solve_expect_order',4,...
    'optimset',optimset('MaxNodes',20,'TolFun',sqrt(eps),...
    'MaxTime',3*60,'MaxIter',2000));
mc3=mc3.estimate('optimizer',1);
check_optimum(mc3)

%% constant parameter model with conditioning information on shorter 
% horizon for the conditioning data and different hypothesis
% add conditional data, change the horizon to 4. In this case all shocks
% are anticipated 4 periods in advance
mc4=mc.set_options('cond_data_ct',cond_data.window('','','',1:4),...
    'solve_expect_order',4,'forecast_conditional_hypothesis','iwb',...
    'optimset',optimset('MaxNodes',20,'TolFun',sqrt(eps),...
    'MaxTime',3*60,'MaxIter',2000));
mc4=mc4.estimate('optimizer','bee_gate');
check_optimum(mc4)

%% constant parameter model with conditioning information in one dimension
% Now we allow only the monetary policy shock (EI) to be anticipated 4
% periods in advance
shock_properties=struct('name','EI','StandardDeviation',nan,'horizon',4);
mc3=mc.set_options('cond_data',cond_data,'shock_properties',shock_properties,...
    'optimset',optimset('MaxNodes',20,'TolFun',sqrt(eps),...
    'MaxTime',3*60,'MaxIter',2000));
mc3=mc3.estimate('optimizer','bee_gate');
check_optimum(mc3)

%% One-Markov-chain model object
profile on
OneChain=rise('Canonical_1_MarkovChain','data',data,...
    'optimset',optimset('MaxNodes',20,'TolFun',sqrt(eps),'MaxTime',15*60,'MaxIter',2000),...
    'check_stability',false,'solver','functional_iteration');%
OneChain=OneChain.estimate('optimizer','bee_gate');
check_optimum(OneChain)
profile off
profile viewer

%% Two-Markov-chain model object
profile on
TwoChain=rise('Canonical_2_MarkovChains','data',data,...
    'optimset',optimset('MaxNodes',20,'TolFun',sqrt(eps),'MaxTime',2*60*60,...
    'MaxIter',2000),'check_stability',false,'solver','functional_iteration');%newton_kronecker_iteration
TwoChain=TwoChain.estimate('optimizer','bee_gate');%
% TwoChain=TwoChain.estimate(1);%
check_optimum(TwoChain)
profile off
profile viewer
%%
close
save_under0='TwoChain_smoothprob';
hfig=figure('name','smoothed probabilities');
subplot(2,1,1)
plot(TwoChain.Filters.smoothed_probabilities.a('a_1'),'linewidth',1.5)
title('chain a, state 1','FontSize',15)
subplot(2,1,2)
plot(TwoChain.Filters.smoothed_probabilities.b('b_1'),'linewidth',1.5)
title('chain b, state 1','FontSize',15)
[~,time_handle]=sup_label('time','x');
[~,prob_handle]=sup_label('probability','y');
set(time_handle,'fontsize',15)
set(prob_handle,'fontsize',15)
% saveas(hfig,[save_under0,'.pdf'])
% saveas(hfig,[save_under0,'.fig'])
% saveas(hfig,[save_under0,'.eps'])

%    [ax1,h1]=sup_label('super X label');
%    [ax2,h2]=sup_label('super Y label','y');
%    [ax3,h2]=sup_label('super Y label (right)','yy');
%    [ax4,h3]=sup_label('super Title'  ,'t');
%    set(h3,'FontSize',30)
%% construct the dynare conference model object
dc=rise('Canonical','data',data,'solver','functional_iteration',...
    'check_stability',false);
dc=dc.estimate('optimizer',1);
check_optimum(dc)
%% close pool of workers
% matlabpool close

