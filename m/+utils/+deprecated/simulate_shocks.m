function [y1,state]=simulate_shocks(y0,T,R,PAI)
% simulate shocks and keeps the effects separated. 

[endo_nbr,exo_nbr,horizon,nstates]=size(R);

if nargin<4
% PAI=Q'*PAI;
    PAI=1/nstates*ones(nstates,1);
end

if isempty(y0)
    y0=zeros(endo_nbr,exo_nbr);
else
    if ~isequal(size(y0),[endo_nbr,exo_nbr])
        error(['y0 expected to be of size [',int2str(endo_nbr),',',int2str(exo_nbr),']'])
    end
end

% draw a state
csp=[0;cumsum(PAI)];
csp(end)=1;
state=find(csp>rand,1,'first')-1;
% for efficiency, exploit sparsity
statecols=any(T(:,:,state));
% unconditional forecasts
y1=T(:,statecols,state)*y0(statecols,:);
% add the shocks
for hh=1:horizon
    shocks=randn(1,exo_nbr);
    y1=y1+R(:,:,hh,state).*shocks(ones(endo_nbr,1),:);
end

