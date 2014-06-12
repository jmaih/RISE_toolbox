function [T,steadyState,state_vars_location]=load_order_var_solution(obj)

% all variables are state variables in the VAR. So we set the location to
% empty
state_vars_location=[];

% we obtain the steady state from the constant of the model


% T is a cell containing in columns, the reduced-form solution of all the
% regimes.
% Each cell contains (in this specific order)
% 1- the autoregressive part
% 2- a column for the risk, which is zero in this case
% 3- the coefficients on the shocks

error('this method is not ready')

% 
% regimes_number=obj.markov_chains.regimes_number;
% order=obj.options.solve_order;
% T=cell(order,regimes_number);
% zzz=repmat('z',1,order);
% ov=obj.order_var.after_solve;
% state_vars_location=obj.locations.after_solve.t.pb;
% steadyState=cell(1,regimes_number);
% for io=1:order
%     for ireg=1:regimes_number
%         if io==1
%             steadyState{ireg}=obj.solution.ss{ireg}(ov);
%         end
%         T{io,ireg}=obj.solution.(['T',zzz(1:io)]){ireg}(ov,:);
%     end
% end

end