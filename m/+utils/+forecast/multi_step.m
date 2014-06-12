function [sims,states,retcode]=multi_step(y0,ss,T,shocks,states,Q,PAI,options)

endo_nbr=size(y0.y,1);

sims=nan(endo_nbr,options.nsteps);

[states(1),Q,PAI,retcode]=generic_tools.choose_state(states(1),Q,PAI,y0.y);

for t=1:options.nsteps
    
    if ~retcode
        
        rt=states(t);
        
        state_vars_location=T{end,rt}; % last row includes the state variables
        
        y1=utils.forecast.one_step(T(:,rt),y0,ss{rt},state_vars_location,...
            options.simul_sig,shocks(:,t+(0:options.k_future)),options.simul_order);
        
        if t<options.nsteps
            [states(t+1),Q,PAI,retcode]=generic_tools.choose_state(states(t+1),Q,PAI,y1.y);
        end
        
        if t>options.burn
            sims(:,t-options.burn)=y1.y;
        end
        
        y0=y1;
    end
    
end

end