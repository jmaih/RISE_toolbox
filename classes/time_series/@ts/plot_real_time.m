function plot_handle=plot_real_time(varargin)
% plot_real_time - hairy plot
%
% Syntax
% -------
% ::
%
%   plot_handle=plot_real_time(rts)
%   plot_handle=plot_real_time(xrange,rts)
%   plot_handle=plot_real_time(rts,varargin)
%   plot_handle=plot_real_time(xrange,rts,varargin)
%
% Inputs
% -------
%
% - **rts** [ts|rts]: valid time series object with possibly several
% columns and exactly one page
%
% - **xrange** []: range over which to restrict the plots
%
% - **varargin** : valid matlab arguments for plot, coming in pairs
%
% Outputs
% --------
%
% - **plot_handle** : handle to the plot
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:


plot_handle0=utils.plot.myplot(@plot_real_time,varargin{:});
if nargout
    plot_handle=plot_handle0;
end

end
