%--- help for generic/draw_parameter ---
%
%  Random parameter draws for RISE model objects.
% 
%  ::
% 
%    [draw,obj] = draw_parameter(obj,simulation_folder);
% 
%  Args:
% 
%     obj (rise | dsge | rfvar | svar): RISE model object
%     simulation_folder (char | struct):
% 
%       - char(1): simulation folder : the stored elements should be structures
%         with fields:
% 
%           - **x** : parameter vectors
%           - **f** : value of minus(log posterior kernel)
% 
%       - char(2): ['mode'|'prior']: draw from the prior distribution or from a
%         multivariate normal distribution around the mode.
% 
%  Returns:
%     :
% 
%     - **draw** [cell]: the first entry is the names of the estimated
%       parameters and the second is a vector of drawn parameters. The whole cell
%       can be pushed in to a model as obj=set(obj,'parameters',draw).
% 
%     - **obj** [rise|dsge|rfvar|svar]: RISE model object in which the drawn
%       parameter has been pushed.
% 
%