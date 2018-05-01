function [sim1,regimes,Qt,retcode]=simul_regime_switch(y0,T,ss,...
state_vars_location,options)

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

simul_sig=options.simul_sig;

simul_order=options.simul_order;

PAI=options.PAI;

Qfunc=options.Qfunc;

endo_nbr=size(y0.y,1);

switch_rule=y0.switch_rule;

inv_order_var=y0.inv_order_var;

shocks=y0.econd.data;

k_future=options.k_future;

% use first page only as it is the most relevant
regimes=vec(y0.rcond.data(:,:,1));

span=options.nsteps+options.burn;

retcode=0;

y0.econd.data=shocks(:,:,ones(3,1));

if ~iscell(switch_rule)
    
    switch_rule={switch_rule};
    
end

% initialize with nans in order to see if/when simulation breaks down
sim1=nan(endo_nbr,span);

y00=y0;

t=0;

while t<span
    
    t=t+1;
    % compute transition matrix and switching probabilities
    %------------------------------------------------------
    [Q,retcode]=Qfunc(y00.y);
    
    if t==1
        
        Q=full(Q);
        
        Qt=Q(:,:,ones(1,span));
        
    else
        
        Qt(:,:,t)=Q;
        
    end
    
    smallrt=regimes(t);
    
    forced_regime=~isnan(smallrt);
    
    if ~forced_regime
        % the state is not known
        if t==1
            % draw from initial distribution
        else
            % draw conditional on yesterday's state
            PAI=Q(regimes(t-1),:);
            
        end
        
    end
    
    shkt=shocks(:,t+(0:k_future));
    
    [y1,regimes(t)]=simulate_forecast(y00,shkt,PAI,smallrt,regimes(t));
    
    sim1(:,t)=y1.y;
    
    % initial conditions for next step
    y00=y1;
    
end

    function [y1,rt]=simulate_forecast(y0,shkt,PAI,rt,reference_rt)
        
        y1=[];
        
        state_list=1:h;
        
        if ~isempty(switch_rule)
            
            [y1,rt]=apply_switch_rule(y0,rt);
            
            return
            
        end
        
        ok=false;
        
        while ~ok
            
            [cp,retcode]=rebuild_cp();
            
            if retcode,return, end
            
            if isnan(rt)
                
                lucky=find(cp>rand,1,'first')-1;
                
                rt=state_list(lucky);
                
            end
            
            myorder={T(:,rt),y0,ss{rt},state_vars_location,...
                simul_sig,shkt,simul_order};
            
            % use the solution prescribed by simul_anticipate_zero
            y1=utils.forecast.one_step_fbs(myorder{:});
            % y1=utils.forecast.one_step_engine(T,y0,ss,xloc,sig,shocks,order)
            
            ok = isempty(options.complementarity)||...
                ~options.complementarity(y1.y);
            
            if ~ok
                
                state_list(state_list==rt)=[];
                                
                if forced_regime
                    
                    warning(['forced to apply solution %0.0f ',...
                        'but the path is infeasible'],reference_rt)
                    
                end
                
                if forced_regime||isempty(state_list)
                    
                    retcode=703;
                    
                    return
                    
                end
                
                rt=nan;
                
            end
            
        end
            
        function [cp,rcode]=rebuild_cp()
            
            rcode=0;
            
            PAI00=PAI(state_list);
            
            PAI00=PAI00/sum(PAI00);
            
            if any(isnan(PAI00))
                
                rcode=703;
                
            end
            
            cp=cumsum(PAI00);
            
            cp=[0,cp(:).'];
            
        end
        
        function [y1,smallrt]=apply_switch_rule(y0,smallrt)
            
            % evaluate all regimes
            %---------------------
            ysr=y0(ones(1,h));
            
            for ireg=1:h
                
                ysr(ireg)=utils.forecast.one_step_fbs(T(:,ireg),...
                    y0,ss{ireg},state_vars_location,...
                    simul_sig,shkt,simul_order);
                % y1=utils.forecast.one_step_engine(T,y0,ss,xloc,sig,shocks,order)
                
            end
            
            ysry=[ysr.y];
            
            good=switch_rule{1}(ysry(inv_order_var,:),smallrt,...
                regimes(1:t-1),...
                sims(inv_order_var,1:t-1),switch_rule{2:end});
            
            % zero probability for offenders
            %-------------------------------
            PAI(~good)=0;
            
            cp=rebuild_cp();
            
            if isnan(smallrt)
                
                lucky=find(cp>rand,1,'first')-1;
                
                smallrt=state_list(lucky);
                
            end
            
            y1=ysr(smallrt);
            
        end
        
    end

end