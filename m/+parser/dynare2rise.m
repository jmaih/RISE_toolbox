%  dynare2rise -- converts a basic dynare file into a RISE one
% 
%  ::
% 
%    duplicates=dynare2rise(dynFileName)
%    duplicates=dynare2rise(dynFileName,riseFileName)
%    duplicates=dynare2rise(dynFileName,riseFileName,stderr_name)
%    duplicates=dynare2rise(dynFileName,riseFileName,stderr_name,detrend)
% 
%  Args:
% 
%     - **dynFileName** [char] : name of the dynare model file with or
%       without extension (.dyn or .mod)
% 
%     - **riseFileName** [char|{'dynFileName.rs'}] : name of the created
%       RISE file 
% 
%     - **stderr_name** [char|{'std'}] : prepended name for the newly
%       created parameters (RISE transforms all the variances and standard
%       deviations into parameters)
% 
%     - **detrend** [true|{false}]: 
%       - when true, the model is written in stationary form like dynare
%         would do behind the scenes. For log growth variables,
%         y-->(y*trend). For linear growth variables, y-->(y+trend). The
%         steady state model in dynare then corresponds exactly to the
%         steady state file produced by RISE, except that the trend
%         variables are added as new variables. 
%       - when false, the model is kept in nonstationary form. For log
%         growth variables the growth rate in an expression of the form f(x)
%         becomes log(f(exp(x)). The growth rate for linear growth variables
%         remains unchanged. The variables detected to have steady state = 0
%         but that need to be stationarized are assigned a growth rate of 0
%         and are declared as LEVEL variables. The steady state remains
%         unchanged.
%       - The steady state for all trend variables with log growth is
%         normalized to 1, while the steady state for trend variables with
%         linear growth is normalized to 0.
% 
%  Returns:
%     :
% 
%     - **duplicates** [struct] : structure containing information on the
%       redundancies
% 
%  NB: the created files are 
% 
%     - riseFileName.rs 
% 
%     - riseFileName_params.m : parameter function file with calibration and
%       priors as separate outputs
% 
%     - riseFileName_sstate.m (optional): steady state function file
% 
%  See also :
%