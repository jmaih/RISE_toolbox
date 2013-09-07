function [x,fval,exitflag,obj]=cmsa_gate_conclude(PROBLEM)

[x,~,~,obj]=cmsa_gate(PROBLEM);

Objective=PROBLEM.objective;
lb=PROBLEM.lb;
ub=PROBLEM.ub;
options=PROBLEM.options;

if exist('fmincon.m','file')
	[x,fval,exitflag]=fmincon(Objective,x,[],[],[],[],lb,ub,[],options);
else
	[x,fval,obj]=local_optimize(Objective,x,fval,lb,ub,options);
	exitflag=0; % disp('Too many function evaluations or iterations.')
end
                    
