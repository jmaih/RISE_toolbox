function [shocks]=create_shocks(nx,shock_id,det_vars,options,user_input)
% INTERNAL FUNCTION
%

if nargin<5
    
    user_input=[];
    
end

% The shock will hit either in the current period,horizon_hit=0 or in a
% later period horizon_hit>0. In principle horizon_hit can differ from the
% number of periods agents can see into the future.

random=options.simul_shock_uncertainty;
% stochastic shocks
%------------------
shocks_span=options.nsteps+options.k_future+options.burn;

shocks=zeros(nx,shocks_span);

if random
    
    set_all_shocks()
    
    shocks=utils.forecast.nullify_deterministic_shocks(shocks,det_vars);
    
end

if ~isempty(shock_id)
    
    impulse=options.impulse;
    
    utils.forecast.check_shock_id(shock_id,nx);
    % first batch/period of shocks
    %-----------------------------
    shocks(shock_id,options.k_future+1)=impulse;
    
end

    function set_all_shocks()
        
        if isempty(user_input)
            
            shocks=randn(nx,shocks_span);
            
            return
            
        end
        
        logic=true(1,nx);
        
        forbidden=cell2mat(user_input(:,1));
        
        logic(forbidden)=false;
        
        shocks(logic,:)=randn(sum(logic),shocks_span);
        
        nuser=numel(forbidden);
        
        for ii=1:nuser
            
            position=forbidden(ii);
            
            [routine,vargs]=utils.code.user_function_to_rise_function(user_input{ii,2});
            
            shocks(position,:)=routine(shocks_span,vargs{:});
            
        end
        
    end

end