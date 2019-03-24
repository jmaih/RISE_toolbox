%--- help for ts/plot_real_time ---
%
%  plot_real_time - hairy plot
% 
%  ::
% 
%     plot_handle=plot_real_time(rts)
%     plot_handle=plot_real_time(xrange,rts)
%     plot_handle=plot_real_time(rts,varargin)
%     plot_handle=plot_real_time(xrange,rts,varargin)
% 
%  Args:
% 
%     rts (ts | rts): valid time series object with possibly several
%       columns and exactly one page
%     xrange : range over which to restrict the plots
%     varargin : valid matlab arguments for plot, coming in pairs
% 
%  Returns:
%     :
% 
%     - **plot_handle** : handle to the plot
% 
%