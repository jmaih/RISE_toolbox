function [x1,x1_linear]=simulation_engine(obj,...
    x0,... % initial conditions
    x0_linear,... % initial conditions for pruning
    z,... % deterministic variables
    regime,...
    shocks)
persistent m_x m_e simul_with_shocks is_observed exo_nbr nlags endo_nbr
if nargin<6
    shocks=[];
end
if isempty(x0_linear)
    nlags=obj.nlags;
    h=obj.markov_chains.regimes_number;
    [m_x,~,m_e]=vartools.resolve(obj.solution,nlags,h);
    simul_with_shocks=~obj.options.simul_no_shocks;
    is_observed=obj.exogenous.is_observed;
    exo_nbr=sum(obj.exogenous.number);
    endo_nbr=obj.endogenous.number(end);
end
if isempty(shocks)
    shocks=generic_tools.set_exogenous_data(exo_nbr,is_observed,simul_with_shocks,z);
end
% autoregressive part
%--------------------
x0_hat=x0(:,end:-1:1);
By=m_x{regime}(:,1:nlags*endo_nbr);
x1=By*x0_hat(:);

% deterministic variables
%------------------------
Bz=m_x{regime}(:,nlags*endo_nbr+1:end);
if obj.constant
    x1=x1+Bz(:,end);
    Bz(:,end)=[];
end
% shocks and deterministic variables
%-----------------------------------
if size(shocks,2)
    Be=zeros(endo_nbr,exo_nbr);
    Be(:,~is_observed)=m_e{regime}(:,:,1);
    Be(:,is_observed)=Bz;
    x1=x1+Be*shocks(:,1);
end
x1=[x0(:,2:end),x1];
x1_linear=x1;
end
