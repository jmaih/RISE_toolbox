%  plots the curvature of an estimated parameter at the model
% 
%  ::
% 
%    [h,legend_]=curvature(xname,db)
% 
%  Args:
% 
%     - **xname** [char]: name of the parameter to plot
%     - **db** [struct]: structure containing the various parameters. Each
%       parameter field is itself a structure with the following
% 
%       - **tex_name** [char]: name of the parameter as it should appear in the
%         title of the plot
%       - **mode** [scalar]: value of the parameter at the mode
%       - **log_post_mode** [scalar]: value of the log posterior mode
%       - **log_lik_mode** [scalar]: value of the log likelihood at the mode
%       - **x** [vector]: x-axis values
%       - **log_post** [vector]: value of the log-posterior for each value of x
%       - **log_lik** [vector]: value of the log-likelihood for each value of x
% 
%  Returns:
%     :
% 
%     - **h** [handle]: handle for the plot
%     - **legend_** [cellstr]: names of the lines in the plot
%     - **tex_name** [char]: name of the parameter
% 
%