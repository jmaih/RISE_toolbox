function contrib=nonlinear_shock_decomp(fsim,smoothed,shocks,nsim,varargin)

y0=smoothed(:,1);

nshocks=size(shocks,1);

contrib=smoothed(:,:,ones(nshocks+1,1));

contrib(:)=0;

% first period = initial conditions only
%----------------------------------------
contrib(:,1,end)=y0;

for ishk=1:nshocks
    
    shocks_i=shocks(:,2:end);
    
    shocks_i(ishk,:)=0;
    
    Si=utils.filtering.nonlinear_counterfactual(y0,fsim,shocks_i,nsim,varargin{:});
        
    contrib(:,2:end,ishk)=smoothed(:,2:end)-Si;
    
end

% initial condition
%--------------------
contrib(:,2:end,end)=smoothed(:,2:end)-sum(contrib(:,2:end,1:end-1),3);

end