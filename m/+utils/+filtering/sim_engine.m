function [s,regimes]=sim_engine(y0,shks,ss,T,Qfunc,isstate,pai0)

% transition matrix based on initial conditions: At this stage The function
% already accounts for the fact that y0 will be re-ordered
%--------------------------------------------------------------------------
Q=Qfunc(y0);

smpl=size(shks,2);

s=y0(:,ones(smpl,1));

% where we were yesterday is probabilistic
%-----------------------------------------
rlag=draw_regime(pai0);

regimes=nan(1,smpl);

order=size(T,1);

for t=1:smpl
    % draw the next regime given the regime before
    %----------------------------------------------
    pai=Q(rlag,:);
    
    rt=draw_regime(pai);
    
    regimes(t)=rt;
    
    % compute the forecast
    %---------------------
    y1=compute_forecast(rt);
    
    s(:,t)=y1;
    
    % initial conditions for next iteration
    %--------------------------------------
    y0=y1;
    
    Q=Qfunc(y0);
    
    rlag=rt;
    
end

    function y=compute_forecast(rt)
        
        y0t=y0-ss{rt};
        
        zt=[y0t(isstate);1;shks(:,t)];
        
        y=ss{rt};
        
        zkron=zt;
        
        for io=1:order
            
            y=y+T{io,rt}*zkron;
            
            if io<order
                
                zkron=kron(zkron,zt);
                
            end
            
        end
        
    end

    function r=draw_regime(pai)
        
        cp=[0,cumsum(pai(:).')];
        
        u=rand;
        
        r=find(cp>u,1,'first')-1;
        
    end

end