function obj=setup_calibration(obj,Calibration)
% setup_calibration -- set parameters
%
% Syntax
% -------
% ::
%
%   obj=setup_calibration(obj,Calibration)
%
% Inputs
% -------
%
% - **obj** [rise_generic]: model object
%
% - **Calibration** [struct|cell]: calibration to push. There are two
% possibilities:
%   - Calibration is a struct: the fields are the parameter names and each
%   name contains a parameter value.
%   - Calibration is a cell: the cell has two columns. The first column
%   holds the names of the parameters to change and the second column their
%   values.
%
% Outputs
% --------
%
% - **obj** [rise|dsge]: model object
%
% More About
% ------------
%
% - the parameter values could be vectors if e.g. we want to take the entire
% parameterization of one model and push it into an identical model. But it
% is not allowed to have situations where one parameter is a vector and
% some other is not. Then RISE will complain that the parameter is not
% controlled by the const markov chain.
%
% Examples
% ---------
%
% See also: 

if isempty(Calibration)
    return
end

obj=setup_calibration@rise_generic(obj,Calibration);
% ensure the model will be re-solved no matter what
%---------------------------------------------------
obj.warrant_resolving = true;

end