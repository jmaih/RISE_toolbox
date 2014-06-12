function ovSolution=load_order_var_solution(obj,y0)

regimes_number=obj.markov_chains.regimes_number;
order=obj.options.solve_order;
T=cell(order+1,regimes_number);
zzz=repmat('z',1,order);
ov=obj.order_var.after_solve;
state_vars_location=obj.locations.after_solve.t.pb;
steady_state=cell(1,regimes_number);
for io=1:order
    for ireg=1:regimes_number
        if io==1
            steady_state{ireg}=obj.solution.ss{ireg}(ov);
            T{end,ireg}=state_vars_location;
            % ov the initial conditions as well
            %----------------------------------
            y0(ireg).y=y0(ireg).y(ov,:);
        end
        T{io,ireg}=obj.solution.(['T',zzz(1:io)]){ireg}(ov,:);
    end
end

ovSolution=struct('T',{T},...
    'y0',{y0},...
    'steady_state',{steady_state});
end