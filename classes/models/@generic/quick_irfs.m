%--- help for quick_irfs ---
%
%  Quick plotting of impulse response functions
% 
%  ::
% 
% 
%    hdl=quick_irfs(m,myirfs)
%    hdl=quick_irfs(m,myirfs,var_list)
%    hdl=quick_irfs(m,myirfs,var_list,shock_list)
%    hdl=quick_irfs(m,myirfs,var_list,shock_list,r0c0)
%    hdl=quick_irfs(m,myirfs,var_list,shock_list,r0c0,xrange)
%    hdl=quick_irfs(m,myirfs,var_list,shock_list,r0c0,xrange,suplab)
% 
%  Args:
% 
%     m (rise | dsge | svar | rfvar): model object
% 
%     myirfs (struct): structure with irfs
% 
%     var_list (cellstr): list of variables of interest
% 
%     shock_list (cellstr): list of shocks of interest
% 
%     r0c0 (vector | {[4,4]}): number of rows and columns in figure
% 
%     xrange (serial | char | {[]}): range over which to plot
% 
%     suplab (false | {true}): add a sup-label to the figure or not
% 
%  Returns:
%     :
% 
%     - **hdl** [handle]: handle to the plotted objects
% 
%  See also: quick_plots
% 
%