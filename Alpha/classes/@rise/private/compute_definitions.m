function obj=compute_definitions(obj,pp)

definitions=obj.func_handles.definitions;
if isempty(obj.solution)
    obj.solution=rise.initialize_solution_or_structure('solution',obj.markov_chains.regimes_number);
end
number_of_regimes=obj.markov_chains.regimes_number;
if nargin<2
    pp=obj.parameter_values;
end
% evaluate definitions
for ii=number_of_regimes:-1:1 % backward to optimize speed
    obj.solution.definitions{ii}=online_function_evaluator(definitions,pp(:,ii)); %#ok<*EVLC>
end
