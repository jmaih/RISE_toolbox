function contrib=nonlinear_shock_decomp(fsim,smoothed,shocks,pos_groups,nsim,varargin)

y0=smoothed(:,1);

nshocks=numel(pos_groups);

contrib=smoothed(:,:,ones(nshocks+1,1));

contrib(:)=0;

% first period = initial conditions only
%----------------------------------------
contrib(:,1,end)=y0;

for ishk=1:nshocks
    
    whichOnes=pos_groups{ishk};
    
    shocks_i=shocks(:,2:end);
    
    shocks_i(whichOnes,:)=0;
    
    Si=utils.filtering.nonlinear_counterfactual(y0,fsim,shocks_i,nsim,varargin{:});
        
    contrib(:,2:end,ishk)=smoothed(:,2:end)-Si;
    
end

% initial condition
%--------------------
contrib(:,2:end,end)=smoothed(:,2:end)-sum(contrib(:,2:end,1:end-1),3);

end