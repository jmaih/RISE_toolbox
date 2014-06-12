function [shocks]=create_shocks(nx,shock_id,det_vars,options)
% function [shocks,State,Q,PAI,retcode]=create_shocks(nx,shock_id,options.k_future,det_vars,random,impulse,...
%     State,Q,PAI,y0)
% The shock will hit either in the current period,horizon_hit=0 or in a later period horizon_hit>0
% horizon_hit can differ from the order at which the model is solved

random=options.random;
impulse=options.impulse;
% stochastic shocks
%------------------
if random
    shocks=randn(nx,options.nsteps+options.k_future);
    shocks=utils.forecast.nullify_deterministic_shocks(shocks,det_vars);
else
    shocks=zeros(nx,options.nsteps+options.k_future);
end

if ~isempty(shock_id)
    utils.forecast.check_shock_id(shock_id,nx);
    % first batch/period of shocks
    %-----------------------------
    shocks(shock_id,options.k_future+1)=impulse;
end

% % markov chains
% %--------------
% [State,Q,PAI,retcode]=generic_tools.choose_state(State,Q,PAI,y0(:,end));

end