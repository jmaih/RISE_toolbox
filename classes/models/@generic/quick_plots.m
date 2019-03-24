%--- help for quick_plots ---
%
%  quick plots
% 
%  ::
% 
%    hdl=quick_plots(m,batch)
%    hdl=quick_plots(m,batch,var_list)
%    hdl=quick_plots(m,batch,var_list,fig_title)
%    hdl=quick_plots(m,batch,var_list,fig_title,r0c0)
%    hdl=quick_plots(m,batch,var_list,fig_title,r0c0,xrange)
% 
%  Args:
% 
%     m (rise | dsge | svar | rfvar): model object
% 
%     batch (struct): structure with variables to plot
% 
%     var_list (cellstr): list of variables of interest
% 
%     fig_title (char | {'no title'}): title of the figure
% 
%     r0c0 (vector | {[4,4]}): number of rows and columns in figure
% 
%     xrange (serial | char | {[]}): range over which to plot
% 
%  Returns:
%     :
% 
%     - **hdl** [handle]: handle to the plotted objects
% 
%  See also:
%     - quick_irfs
% 
%