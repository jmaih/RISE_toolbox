%--- help for dsge/fold_solution ---
%
%  fold_solution : folds higher-order solutions in order to facilitate
%  comparison with other software
% 
%   Syntax :
% 
%    [foldsol,zlist_reordered,ylist_reordered]=fold_solution(m)
% 
%    [foldsol,zlist_reordered,ylist_reordered]=fold_solution(m,order)
% 
%    [foldsol,zlist_reordered,ylist_reordered]=fold_solution(m,order,user_state_list)
% 
%    [foldsol,zlist_reordered,ylist_reordered]=fold_solution(m,order,user_state_list,user_endo_list)
% 
%   Inputs:
% 
%   - **m** [rise|dsge]: scalar or vector of solved model objects, possibly
%     with multiple parameterizations
% 
%   - **order** [numeric|empty]: Desired order of perturbation <= order of solution
% 
%   - **user_state_list** [cellstr|empty]: List of state variables as
%     appearing in the non-RISE solution
% 
%   - **user_endo_list** [cellstr|empty]: List of endogenous variables in the
%     non-RISE solution
% 
%   Outputs:
% 
%   - **foldsol** [cell array]: folded solutions
% 
%   - **zlist_reordered** [cellstr]: reordered list of state variables
% 
%   - **ylist_reordered** [cellstr]: reordered list of endogenous variables
%