%--- help for mcmc/statistics ---
%
%  statistics : statistics from the mcmc chains
% 
%  ::
% 
%     s=statistics(obj);
% 
%     s=statistics(obj,opts);
% 
%  Args:
% 
%     - obj (mcmc object): mcmc object
% 
%     - opts (struct): structure with options of which statistics to compute.
%       for more information see help utils.mcmc.statistics
% 
%  Returns:
%     :
% 
%     - **s** (struct): structure with statistics of interest
% 
%     - **statTab** (table): table representation of some of the statistics
%