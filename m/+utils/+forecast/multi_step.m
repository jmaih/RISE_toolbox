function [sims,states,retcode,Qt]=multi_step(y0,ss,T,state_vars_location,options)
% Qt below may be used to get the time series of the transition matrix
endo_nbr=size(y0.y,1);
PAI=options.PAI;
Qfunc=options.Qfunc;
states=options.states;
shocks=options.shocks;
options=rmfield(options,{'states','shocks','PAI','Qfunc','y'});

sims=nan(endo_nbr,options.nsteps);

[states(1),Q0,PAI,retcode]=generic_tools.choose_state(states(1),Qfunc,PAI,y0.y);
do_Qt=nargout>3;
span=options.nsteps+options.burn;
if do_Qt
    Qt=Q0(:,:,ones(1,options.nsteps));
end
for t=1:span
    
    if ~retcode
        if t==options.burn && do_Qt
            Qt(:,:,t-options.burn+1)=Q0;
        end
        rt=states(t);
        
        y1=utils.forecast.one_step(T(:,rt),y0,ss{rt},state_vars_location,...
            options.simul_sig,shocks(:,t+(0:options.k_future)),options.simul_order);
        
        if t<span
            [states(t+1),Q0,PAI,retcode]=generic_tools.choose_state(states(t+1),Qfunc,PAI,y1.y);
        end
        
        if t>options.burn
            sims(:,t-options.burn)=y1.y;
            if do_Qt && t<span
                Qt(:,:,t-options.burn+1)=Q0;
            end
        end
        y0=y1;
    end
end

states=states(options.burn+1:end);

end