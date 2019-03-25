%  linear_dsge_generator -- generates a linear constant-parameter dsge model
% 
%  ::
% 
%    linear_dsge_generator(riseFile)
%    linear_dsge_generator(riseFile,endo_list)
% 
%  Args:
% 
%     riseFile (char): name of the rise file to be generated
% 
%     endo_list (cellstr | char| numeric|{3}): list or number of endogenous
%       variables 
% 
%  Returns:
%     :
% 
%     nothing
% 
%  Note:
% 
%     - the endogenous variables are denoted by "x" if the list of
%       endogenous variables is not provided
% 
%     - the exogenous variables are denoted by "e"
% 
%     - the coefficients on leads of endogenous variables by "ap_i_j"
% 
%     - the coefficients on current endogenous variables by "a0_i_j"
% 
%     - the coefficients on lags of endogenous variables by "am_i_j"
% 
%     - the constant terms by "b_i"
% 
%     - the coefficients on exogenous variables by "c_i"
% 
%     - where "i" denotes the equation, "j" the variable
% 
%