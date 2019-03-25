%  plots prior and, if information is available, posterior densities
% 
%  ::
% 
%    [h,legend_]=prior_posterior(ss,varargin)
% 
%  Args:
% 
%     - **ss** [struct]: structure containing the relevant information to plot
%       and more specifically
% 
%       - **x_kdens** [vector]: x-axis values for the posterior density
%       - **f_kdens** [vector]: y-axis values for the posterior density
%       - **x_prior** [vector]: x-axis values for the prior density
%       - **f_prior** [vector]: y-axis values for the prior density
%       - **mean_sim** [scalar]: mean of the posterior simulation
%       - **post_mode** [scalar]: value at the posterior maximization mode
%       - **post_mode_sim** [scalar]: value at the posterior simulation mode
%       - **tex_name** [char]: name of the parameter
% 
%     - **varargin** [pairwise arguments]: standard plotting arguments for
%       matlab
% 
%  Returns:
%     :
% 
%     - **h** [handle]: handle for the plot
%     - **legend_** [cellstr]: names of the lines in the plot
%     - **tex_name** [char]: name of the parameter
% 
%  Note:
% 
%     - Only the prior density is plotted if no posterior information is
%       available.
% 
%