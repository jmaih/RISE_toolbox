function list=input_list()
% INTERNAL FUNCTION
%

list={
    'y',... endogenous
    'x',... exogenous
    'ss',... steady state
    'param',... parameters
    'def',... definitions
    's0',... current regime
    's1' ... next period's regime
    };
end