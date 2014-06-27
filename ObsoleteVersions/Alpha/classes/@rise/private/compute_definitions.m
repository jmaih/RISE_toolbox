function [obj,retcode]=compute_definitions(obj,pp)

definitions=obj.func_handles.definitions;
if isempty(obj.solution)
    obj.solution=rise.initialize_solution_or_structure('solution',obj.markov_chains.regimes_number);
end
number_of_regimes=obj.markov_chains.regimes_number;
if nargin<2
    pp=obj.parameter_values;
end
% evaluate definitions
retcode=0;
for ii=number_of_regimes:-1:1 % backward to optimize speed
    tmp=online_function_evaluator(definitions,pp(:,ii)); 
    if any(isnan(tmp))||any(isinf(tmp))||any(~isreal(tmp))
        retcode=5;
        break
    end
    obj.solution.definitions{ii}=tmp; %#ok<*EVLC>
end
