function h=boxplot(varargin)
% Make box plots of multiple time series in a frame. Overloads
% Matlab's boxplot function for ts objects
%
% ::
%
%    h = boxplot(db);
%    h = boxplot(xrange, db);
%    h = boxplot(..., varargin);
%
% Args:
%
%    - **db** [ts]: time series object
%    - **xrange** [char|cellstr|serial date]: Range of the data to plot
%    - varargin: additional matlab (See BOXPLOT) and RISE (PARSE_PLOT_ARGS)
%      options coming in pairs
%
% See also:
%    parse_plot_args

tmp=utils.plot.myplot(@boxplot,varargin{:});

if nargout

    h=tmp;

end

end
