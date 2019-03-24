%--- help for rfvar/identification ---
%
%  Parses identifying restrictions into RISE interpretable identification function
% 
%  Note:
%     This function is a part of a general workflow. See XXXXXXXXXX to check the workflow.
% 
%  ::
% 
%     Rfunc = identification(var, restr, shock_names, varargin);
% 
%  Args:
%     var (rfvar object): rfvar object
%     restr : restrictions can take two forms
% 
%        - **'choleski'**: choleski identification. Supply the ordering as a cell array as varargin{1}
%        - **cellstr**: n x 2 cell array, with n the number of restrictions.
%           A typical restriction takes the form:
%           - 'vbl@shk','sign'
%           - 'vbl{lags}@shk','sign'
%           where sign = +|-
%           where shk is the name of the shock
%           where vbl is the name of the endogenous variable affected by
%           shock
%           lags with defaut 0 is the list of lags at which the shock
%           affects the variable with a 'sign' sign. lags is any valid
%           expression that matlab can evaluate such as 0, 0:3, [1,4,5],
%           [1,inf,2,5] or more generally a:b, [a:b,c,d], etc.
% 
%     shock_names (cellstr): cell of structural shock names
% 
%  Returns:
%     :
% 
%     - **Rfunc** (function handle): a function that takes parameters and then returns structural shock matrix. This function function is not intended to be used used on its own, and supposed to be used as an argument in different functions.
% 
%  See also:
%     - :func:`autocov <var.autocov>`
%     - :func:`forecast <var.forecast>`
%     - :func:`historical_decomposition <var.historical_decomposition>`
%     - :func:`irf <var.irf>`
%     - :func:`variance_decomposition <var.variance_decomposition>`
% 
%