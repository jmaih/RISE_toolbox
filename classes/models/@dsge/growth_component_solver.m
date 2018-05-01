function T=growth_component_solver(obj,pos,T)
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

is_log_var=obj.endogenous.is_log_var;

ov=obj.order_var;

solve_order=obj.options.solve_order;

zzz=repmat('z',1,solve_order);

myT=cell(solve_order,1);

xloc=pos.z.pb;

nstate_vars=numel(xloc);

k=zeros(obj.endogenous.number,1);

h=numel(T.Tz);

nshocks=size(T.Tz{1},2)-(nstate_vars+1);

zero_shocks=zeros(nshocks,1);

x0=struct('y',[]);

for ireg=1:h
    
    do_one_regime();
    
end

% now push final result back into T: write is as imaginary so as to
% 1- make it readily available if it is to be looked at
% 2- avoid nonlinear computations trying to mix it with Tsig
% 3- give the ability to turn it off separately
% The costs are as follows
% a- hiding it into Tz_sig implies readapting the routines calling
% one_step_engine
% b- how does it square with pruning? A pruned model will not prevent this
% from exploding...
%----------------------------------------------------------------------

    function do_one_regime()
        
        bgp=obj.solution.bgp{ireg};
        
        bgp(is_log_var)=log(bgp(is_log_var));
        
        if ~any(bgp)
            
            return
            
        end
        
        for io=1:solve_order
            
            % No re-ordering of rows at this stage: we are still solving
            %-----------------------------------------------------------
            myT{io}=T.(['T',zzz(1:io)]){ireg};
            
        end
        
        ss=obj.solution.ss{ireg};
        
        ss(is_log_var)=log(ss(is_log_var));
        
        ss=ss(ov);
        
        bgp=bgp(ov);
        
        x0.y=ss;
        
        x1=x0.y+bgp;
        
        sig=0;
        
        x1_k=utils.forecast.one_step_engine(myT,x0,ss,xloc,sig,...
            zero_shocks,solve_order);
        
        k=x1-x1_k.y;
        
        T.Tz{ireg}(:,nstate_vars+1)=T.Tz{ireg}(:,nstate_vars+1)+k*1i;
        
    end

end