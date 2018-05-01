function [shocks]=create_shocks(nx,shock_id,det_vars,options)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:


% The shock will hit either in the current period,horizon_hit=0 or in a
% later period horizon_hit>0. In principle horizon_hit can differ from the
% number of periods agents can see into the future. 

random=options.simul_shock_uncertainty;
% stochastic shocks
%------------------
shocks_span=options.nsteps+options.k_future+options.burn;
if random
    shocks=randn(nx,shocks_span);
    shocks=utils.forecast.nullify_deterministic_shocks(shocks,det_vars);
else
    shocks=zeros(nx,shocks_span);
end

if ~isempty(shock_id)
    impulse=options.impulse;
    utils.forecast.check_shock_id(shock_id,nx);
    % first batch/period of shocks
    %-----------------------------
    shocks(shock_id,options.k_future+1)=impulse;
end

end