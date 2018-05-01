function obj=setup_calibration(obj,Calibration)
% setup_calibration -- set parameters
%
% ::
%
%
%   obj=setup_calibration(obj,Calibration)
%
% Args:
%
%    - **obj** [generic]: model object
%
%    - **Calibration** [struct|cell]: calibration to push. There are two
%    possibilities:
%      - Calibration is a struct: the fields are the parameter names and each
%      name contains a parameter value.
%      - Calibration is a cell: the cell has two columns. The first column
%      holds the names of the parameters to change and the second column their
%      values.
%
% Returns:
%    :
%
%    - **obj** [rise|dsge]: model object
%
% Note:
%
%    - the parameter values could be vectors if e.g. we want to take the entire
%    parameterization of one model and push it into an identical model. But it
%    is not allowed to have situations where one parameter is a vector and
%    some other is not. Then RISE will complain that the parameter is not
%    controlled by the const markov chain.
%
% Example:
%
%    See also:

if isempty(Calibration)
    
    return
    
end

obj=setup_calibration@generic(obj,Calibration);
% ensure the model will be re-solved no matter what
%---------------------------------------------------
obj.warrant_resolving = true;

end