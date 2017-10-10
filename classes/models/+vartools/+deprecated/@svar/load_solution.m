function [T,R,steady_state,new_order,state_vars_location,log_vars]=...
    load_solution(obj,flag,~)
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
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

% when the type is ov, R is empty while T is of size solve_order x
% regimes_number. Basically, T is returned as Tz, Tzz, Tzzz, etc. where the
% rows that were originally in alphabetical order are put in order of
% simulation/processing
% when the type is ov, R corresponds to the shocks and T corresponds to the
% autoregressive elements and the matrices in T are square

if isempty(obj)
    if nargout>1
        error([mfilename,':: number of output arguments cannot exceed 1 when the object is empty'])
    end
    T=struct();
    return
end

log_vars=[];

[T,R,steady_state]=set_solution_to_companion(obj);

nrows=size(T{1},1);
new_order=1:nrows;

state_vars_location=1:nrows;

if strcmp(flag,'ov')
    % we artificially introduce a zero column here in order to be able to
    % use utils.forecast.one_step during simulations and forecasts
    zero_col=zeros(nrows,1);
    for ireg=1:numel(T)
        T{ireg}=[T{ireg},zero_col,R{ireg}];
        R{ireg}=[];
    end
end

end