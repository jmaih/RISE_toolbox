%  aim2rise -- converts a basic dynare file into a RISE one
% 
%  ::
% 
%    duplicates=aim2rise(aimFileName)
%    duplicates=aim2rise(aimFileName,riseFileName)
%    duplicates=aim2rise(aimFileName,riseFileName,stderr_name)
% 
%  Args:
% 
%     - **aimFileName** [char] : name of the dynare model file with or
%       without extension (.aim)
% 
%     - **riseFileName** [char|{'aimFileName.rs'}] : name of the created
%       RISE file 
% 
%     - **stderr_name** [char|{'std'}] : prepended name for the newly
%       created parameters (RISE transforms all the variances and standard
%       deviations into parameters)
% 
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
%  See also : :ref:`dynare2rise`
%