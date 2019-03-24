%--- help for generic/plot_posteriors ---
%
%  Computes posterior densities for estimated parameters
% 
%  ::
% 
%    pdata=plot_posteriors(obj)
% 
%    pdata=plot_posteriors(obj,simfold)
% 
%    pdata=plot_posteriors(obj,simfold,parlist)
% 
%    pdata=plot_posteriors(obj,simfold,parlist,npoints)
% 
%    pdata=plot_posteriors(obj,simfold,parlist,npoints,subset)
% 
%    pdata=plot_posteriors(obj,simfold,parlist,npoints,subset,varargin)
% 
%  Args:
% 
%     obj (rise | dsge | rfvar | svar): model object
% 
%     simfold (empty | char | struct): location of the simulations. If
%       empty, it is assumed that the simulations are saved to disc and are
%       located in the address found in obj.folders_paths.simulations. If it is a
%       "char", this corresponds to the location of the simulation. Otherwise, if
%       it is a struct, then it has to be the output of posterior_simulator.m
% 
%     npoints (numeric | {20^2}): the number of points in the
%       discretization of the prior support
% 
%     parlist (empty | char | cellstr): list of the parameters for which one
%       wants to plot the posteriors
% 
%     subset (cell array|{empty}): When not empty, subset is a
%      1 x 2 cell array in which the first cell contains a vector selecting
%      the columns to retain in each chain and the second column contains
%      the chains retained. Any or both of those cell array containts can be
%      empty. Whenever an entry is empty, all the information available is
%      selected. E.g. subsetting with dropping and trimming
%      mysubs={a:b:c,[1,3,5]}. In this example, the first
%      element selected is the one in position "a" and
%      thereafter every "b" element is selected until we reach
%      element in position "c". At the same time, we select
%      markov chains 1,3 and 5.
% 
%  Returns:
%     :
% 
%     - **pdata** [struct]: optional output argument, pdata is a structure
%       containing the information needed to plot the posterior and prior
%       densities. The user can always plot those using
%       utils.plot.prior_posterior(pdata.(pname)), where pname is the name of
%       one particular parameter of interest.
% 
%     - **hdl** [vector]: optional output argument, placeholder for handles to
%       graphs
% 
%  Note:
% 
%     - if there are no output arguments, figures with posterior and prior
%       marginal densities are plotted, but not saved!!!.
%       see also utils.plot.prior_posterior
% 
%