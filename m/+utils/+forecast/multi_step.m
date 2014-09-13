function [sims,states,retcode]=multi_step(y0,ss,T,state_vars_location,options)
endo_nbr=size(y0.y,1);
PAI=options.PAI;
Q=options.Q;
states=options.states;
shocks=options.shocks;
options=rmfield(options,{'states','shocks','PAI','Q','y'});

sims=nan(endo_nbr,options.nsteps);

[states(1),Q,PAI,retcode]=generic_tools.choose_state(states(1),Q,PAI,y0.y);

span=options.nsteps+options.burn;
for t=1:span
    
    if ~retcode
        
        rt=states(t);
        
        y1=utils.forecast.one_step(T(:,rt),y0,ss{rt},state_vars_location,...
            options.simul_sig,shocks(:,t+(0:options.k_future)),options.simul_order);
        
        if t<span
            [states(t+1),Q,PAI,retcode]=generic_tools.choose_state(states(t+1),Q,PAI,y1.y);
        end
        
        if t>options.burn
            sims(:,t-options.burn)=y1.y;
        end
        
        y0=y1;
    end
    
end

states=states(options.burn+1:end);

end