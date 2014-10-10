function universal_options = optimization_universal_options()
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

% optimization_universal_options: sets the common default options for all
% optimizers.
%   No Detailed explanation needed

universal_options=struct('MaxNodes',20,... % : number of initial solutions [51]
    'MaxIter',4000,...
    'MaxFunEvals',150000,...
    'MaxTime',3600,...
    'start_time',clock,...
    'verbose',10,...
    'stopping_created',false,...
    'nonlcon',[],'TolX',1e-4,'TolFun',1e-4);

end

