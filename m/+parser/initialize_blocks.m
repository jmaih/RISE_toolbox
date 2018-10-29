function [blocks,block_col_names]=initialize_blocks()
% INTERNAL FUNCTION
%

listing=parser.listing();

blocks={
    'endogenous'            ,{'var','endogenous'}         ,listing   ,false
    'exogenous'             ,{'varexo','exogenous'}       ,listing   ,false
    'parameters'            ,'parameters'                 ,listing   ,false
    'observables'           ,{'varobs','observables'}     ,listing   ,false
    'log_vars'              ,{'log_vars','log_variables'} ,listing   ,false
    'level_variables'       ,'level_variables'            ,listing   ,false
    'model'                 ,'model'                      ,cell(0,3) ,false
    'steady_state_model'    ,'steady_state_model'         ,cell(0,3) ,false
    'parameterization'      ,'parameterization'           ,cell(0,3) ,false
    'planner_objective'     ,'planner_objective'          ,cell(0,3) ,false
    'exogenous_definition'  ,'exogenous_definition'       ,cell(0,3) ,false
    'parameter_restrictions','parameter_restrictions'     ,cell(0,3) ,false
    'Legend'                ,'Legend'                     ,[]        ,false
    };

block_col_names={'name','trigger','listing','is_all_but'};

end