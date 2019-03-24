%--- help for generic/mode_curvature ---
%
%  Checks the curvature at the posterior mode
% 
%  ::
% 
%    db = mode_curvature(obj)
% 
%    db = mode_curvature(obj,varlist)
% 
%    db = mode_curvature(obj,varlist,N)
% 
%    db = mode_curvature(obj,varlist,N,type)
% 
%  Args:
% 
%     obj (rise | dsge | rfvar | svar): model object
% 
%     varlist (char | cellstr | empty): list of parameters for which we want
%       to check the curvature
% 
%     N ({20} | integer): Number of grid points
% 
%     type ({'max'} | 'min' | 'range'): normalization of the log-posterior
%       and the log-likelihood.
% 
%     dbin (struct|empty): structure containing the information to plot the
%       curvature. Each field is the name of a particular parameter. This is
%       to avoid a costly recomputation of db
% 
%  Returns:
%     :
% 
%     - **db** [struct|cell array|vector]: structure containing the  
%       information to plot the curvature. Each field is the name of a 
%       particular parameter. Alternatively, when dbin is not empty, db is a 
%       handle to the plots.
% 
%  Note:
% 
%     - when no output is requested, plots are made but not saved.
% 
%     - one way to plot the curvatures from the output is to use the function
%       utils.plot.curvature
% 
%  See also:
%     - utils.plot.curvature
% 
%