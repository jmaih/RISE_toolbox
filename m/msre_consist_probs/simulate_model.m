function [sims,estim_probs,regimes,retcode]=simulate_model(m,inner_probs,x0,...
    T0,shocks,control_shocks,s0,options)

% the user enters x0 in the order of m.endogenous.name

[theSolver,vargs]=utils.code.user_function_to_rise_function(...
    options.msre_solver);


nsteps=size(shocks,3);

sims=nan(numel(x0),nsteps+1);

order_var=m.order_var;

sims(:,1)=x0;

inv_order_var=m.inv_order_var;

regimes=nan(1,nsteps+1);

regimes(1)=s0;

estim_probs=[];

for istep=2:nsteps+1
    
    [T0,retcode,sim_prob,Q]=theSolver(m,inner_probs,T0,sims(:,istep-1),control_shocks,options,vargs{:});
    
    if retcode
        
        return
        
    end
    
    if istep==2
        
        Tsim=T0;
        
    end

    if ~isempty(sim_prob)
        
        if istep==2
            
            nprobs=numel(sim_prob);
            
            estim_probs=nan(nprobs,nsteps);
            
        end
        
        estim_probs(:,istep-1)=sim_prob(:);  %#ok<AGROW>
        
    end
    
    % pick a regime based on Q or something
    %---------------------------------------
    
    probs=[0,cumsum(Q(s0,:))];
    
    probs(end)=1;
    
    s1=find(probs>=rand,1,'first')-1;
    
    regimes(istep)=s1;
    
    Tsim.Tx=T0.Tx(:,:,s1); Tsim.Tsig=T0.Tsig(:,:,s1); Tsim.Te=T0.Te(:,:,s1);
    Tsim.ss=T0.ss(s1);
    
    % change the order here
    %-----------------------    
    tmp=one_step_simulator(Tsim,sims(order_var,istep-1),shocks(:,:,istep-1));
    
    % change again when it comes to storing
    %---------------------------------------
    sims(:,istep)=tmp{1}(inv_order_var);
    
    s0=s1;
    
    fprintf('step %0.0f of %0.0f complete\n\n',istep-1,nsteps);
    
end

