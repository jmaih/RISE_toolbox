%--- help for generic/plot_priors ---
%
%  Computes prior densities for estimated parameters
% 
%  ::
% 
%    ppdata=plot_priors(obj)
% 
%    ppdata=plot_priors(obj,parlist)
% 
%    ppdata=plot_priors(obj,parlist,trunc)
% 
%    ppdata=plot_priors(obj,parlist,trunc,npoints)
% 
%    ppdata=plot_priors(obj,parlist,trunc,npoints,varargin)
% 
%  Args:
% 
%     obj (rise | dsge | rfvar | svar): model object
% 
%     parlist ([] | char | cellstr): list of the parameters for which one
%       wants to plot the priors
% 
%     trunc (numeric | {1e-3}): serves to truncate the support
% 
%     npoints (numeric | {20^2}): the number of points in the
%       discretization of the prior support
% 
%     varargin (pairwise plotting arguments):
% 
%  Returns:
%     :
% 
%     - **pdata** [struct]: optional output argument, pdata is a structure
%       containing the information needed to plot the prior densities. The user
%       can always plot those using utils.plot.prior_posterior(ppdata.(pname)),
%       where pname is the name of one particular parameter of interest.
% 
%     - **hdl** [vector]: optional output argument, placeholder for handles to
%       graphs
% 
%  Note:
% 
%     - if there are no output arguments, figures with prior densities are
%       plotted, but not saved!!!.
% 
%     - There are other arguments that can be set via the model object directly
%       and that are relevant for this function. They are:
% 
%  See also:
%     - utils.plot.plot_posteriors
%     - utils.plot.plot_priors_and_posteriors
% 
%