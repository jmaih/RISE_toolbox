function [sims,states,retcode,Qt]=multi_step(y0,ss,T,state_vars_location,options)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

% Qt below may be used to get the time series of the transition matrix
endo_nbr=size(y0.y,1);
PAI=options.PAI;
Qfunc=options.Qfunc;
states=options.states;
shocks=options.shocks;
simul_update_shocks_handle=options.simul_update_shocks_handle;
simul_do_update_shocks=options.simul_do_update_shocks;
options=rmfield(options,{'states','shocks','PAI','Qfunc','y'});

sims=nan(endo_nbr,options.nsteps);
h=size(T,2);

do_Qt=nargout>3;
span=options.nsteps+options.burn;

retcode=0;
y1=[];
for t=1:span
    if ~retcode
        % process shocks
        %----------------
        shocks_t=shocks(:,t+(0:options.k_future));
        if ~isempty(simul_update_shocks_handle) && simul_do_update_shocks
            shocks_t=simul_update_shocks_handle(shocks_t,y0.y);
        end
        
        % compute transition matrix and switching probabilities
        %------------------------------------------------------
        [Q,retcode]=Qfunc(y0.y);
        if ~retcode
            rt=states(t);
            if isnan(rt)
                if do_Qt
                    if t==1
                        Qt=Q(:,:,ones(1,options.nsteps));
                    end
                    if t>options.burn && t<span
                        Qt(:,:,t-options.burn)=Q;
                    end
                end
                PAI=Q'*PAI;
            else
                % the state is already known
            end
            
            % compute the forecast
            %---------------------
            state_list=1:h;
            iter=0;
            while isempty(y1) && iter < h
                % draw state and compute forecast
                %--------------------------------
                one_step();
                iter=iter+1;
            end
            if isempty(y1)
                error('I could not find a feasible path')
            end
            if isnan(states(t))
                states(t)=rt;
            end
            
            if t>options.burn
                sims(:,t-options.burn)=y1.y;
            end
            y0=y1;
            y1=[];
        end
    end
end
states=states(options.burn+1:end);

    function one_step()
        cp=rebuild_cp();
        if isnan(rt)
            lucky=find(cp>rand,1,'first')-1;
            rt=state_list(lucky);
        end
        y1=utils.forecast.one_step(T(:,rt),y0,ss{rt},state_vars_location,...
            options.simul_sig,shocks_t,options.simul_order);
        if ~options.complementarity(y1.y)
            if ~isnan(states(t))
                error(sprintf('forced to apply solution %0.0f but cannot find a feasible path',states(t))) %#ok<SPERR>
            end
            state_list(state_list==rt)=[];
            rt=nan;
        end
        function cp=rebuild_cp()
            PAI00=PAI(state_list);
            PAI00=PAI00/sum(PAI00);
            if any(isnan(PAI00))
                error('I could not find a feasible path')
            end
            cp=cumsum(PAI00);
            cp=[0,cp(:).'];
        end
    end

end

% function [sims,states,retcode,Qt]=multi_step(y0,ss,T,state_vars_location,options)
% % Qt below may be used to get the time series of the transition matrix
% endo_nbr=size(y0.y,1);
% PAI=options.PAI;
% Qfunc=options.Qfunc;
% states=options.states;
% shocks=options.shocks;
% simul_update_shocks_handle=options.simul_update_shocks_handle;
% simul_do_update_shocks=options.simul_do_update_shocks;
% options=rmfield(options,{'states','shocks','PAI','Qfunc','y'});
% 
% sims=nan(endo_nbr,options.nsteps);
% 
% [states(1),Q0,PAI,retcode]=generic_tools.choose_state(states(1),Qfunc,PAI,y0.y);
% do_Qt=nargout>3;
% span=options.nsteps+options.burn;
% if do_Qt
%     Qt=Q0(:,:,ones(1,options.nsteps));
% end
% for t=1:span
%     
%     if ~retcode
%         if t==options.burn && do_Qt
%             Qt(:,:,t-options.burn+1)=Q0;
%         end
%         rt=states(t);
%         
%         shocks_t=shocks(:,t+(0:options.k_future));
%         if ~isempty(simul_update_shocks_handle) && simul_do_update_shocks
%             shocks_t=simul_update_shocks_handle(shocks_t,y0.y);
%         end
%         
%         y1=utils.forecast.one_step(T(:,rt),y0,ss{rt},state_vars_location,...
%             options.simul_sig,shocks_t,options.simul_order);
%         
%         if t<span
%             [states(t+1),Q0,PAI,retcode]=generic_tools.choose_state(states(t+1),Qfunc,PAI,y1.y);
%         end
%         
%         if t>options.burn
%             sims(:,t-options.burn)=y1.y;
%             if do_Qt && t<span
%                 Qt(:,:,t-options.burn+1)=Q0;
%             end
%         end
%         y0=y1;
%     end
% end
% 
% states=states(options.burn+1:end);
% 
% end