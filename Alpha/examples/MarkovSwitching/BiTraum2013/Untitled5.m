symbols={
'exclamation','!'
'quotation mark','"'
'pound sign','#'
'dollar sign','$'
'percent sign','%'
'and','&'
'apostrophe',''''
'lparenthesis','('
'rparenthesis',')'
'multiplication','*'
'addition','+'
'comma',','
'minus','-'
'rdivide','/'
'integer','0123456789'
'colon',':'
'semicolon',';'
'less','<'
'equal','='
'greater','>'
'question','?'
'at','@'
'alpha','ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
'rbracket','['
'ldivide','\'
'rbracket',']'
'hat','^'
'dash','_'
'lcbrace','{'
'or','|'
'rcbrace','}'
'tilde','~'
};

%%
block_names={
'var',1
'endogenous',1
'varexo',2
'exogenous',2
'shocks',2
'parameters',3
'log_vars',4
'varobs',5
'observables',5
'model',6
'parameterization',7
'steady_state_model',8
'parameter_restrictions',9
'exogenous_definitions',10
'planner_objective',11
};

word_grammar
equation_grammar
model_grammar(equations, definitions, time_varying_probs)
steady_state_model_grammar(equations)
parameterization_grammar
parameter_restrictions_grammar
declarations_grammar
planner_objective_grammar

% +parser classes < handle
% model
% equation
% steady state
% parameterization
% parameter_restrictions
% planner_objective
% declarations

% methods
% grammar
% append

