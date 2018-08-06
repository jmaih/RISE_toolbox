function varargout=plot(varargin)
% Plots a rise time series
%
% ::
%
%    plot(db);
%    plot(daterange,db);
%    plot(daterange,db, varargin);
%
% Args:
%
%    db (ts object): time series object
%    daterange: date range. See ts for different formats
%    varargin: need to come in pairs (except for string format)
%
%       - s (string): line type, symbol spec, and color spec (same as MATLAB's plot function)
%       - 'nticks': integer (default 8) number of xticks
%       - 'date_format': date format option used in MATLAB.
%       - **vline** [char|cellstr|serial dates|{''}] : vertical line(s) e.g.
%         'vline' = '2000Q1'= '2000Q1,2003Q2' must be in the same frequency as
%         the database to be plotted
%       - **hline** [integer|{''}] : horizontal line(s) 'hline' =1, =[1 5.5 2]
%       - **logy** [true|{false}] : log the database or not
%
% Returns:
%    :
%
%      - **varargout** [scalar|vector] : handle to the lines of plot
%
% Example:
%    ::
%
%       plot('1994m7:1997m1',db(:,:,1),...
%           'figsize',[2,2],...
%           'figtitle','no title',...
%           'nticks',10,...
%           'legend',{'v1','v2'},...
%           'legend_loc','BO',...
%           'logy',true,...
%           'secondary_y',{'v1','v4'},...
%           'subplots',true,...
%           'linewidth',2,...
%           'date_format',17);
%       xrotate(90)
%
% Note:
%     In addition to those matlab properties, RISE adds further properties,
%     which allow to control for. See parse_plot_args
%

[varargout{1:nargout}]=utils.plot.myplot(@plot,varargin{:});

end