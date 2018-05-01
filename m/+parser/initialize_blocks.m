function blocks=initialize_blocks()
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

listing=parser.listing();

blocks={
    'endogenous'            ,{'var','endogenous'}         ,listing
    'exogenous'             ,{'varexo','exogenous'}       ,listing
    'parameters'            ,'parameters'                 ,listing
    'observables'           ,{'varobs','observables'}     ,listing
    'log_vars'              ,{'log_vars','log_variables'} ,listing
    'level_variables'       ,'level_variables'            ,listing
    'model'                 ,'model'                      ,cell(0,3)
    'steady_state_model'    ,'steady_state_model'         ,cell(0,3)
    'parameterization'      ,'parameterization'           ,cell(0,3)
    'planner_objective'     ,'planner_objective'          ,cell(0,3)
    'exogenous_definition'  ,'exogenous_definition'       ,cell(0,3)
    'parameter_restrictions','parameter_restrictions'     ,cell(0,3)
    'Legend'                ,'Legend'                     ,[]
    };

end