function plot_handle=barh(varargin)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:


%     barh(X,db) draws the columns of the M-by-N matrix db as M groups of N
%     vertical bars versus vector of dates X. X can be of
%     multiple forms and need not be of the same length as db:
%         - a string of dates. e.g. X='1990:2000' X='1990Q1:2000Q1'
%         X='1990m1:2000m2' 
%         - serial dates. e.g. X=date2serial('1990Q1):date2serial('2000Q1')
%         X='1990m1:2000m2' 
%  
%     barh(db) uses the dates within db.  The colors are set by the colormap.
%  
%     barh(X,db,WIDTH) or barh(db,WIDTH) specifies the width of the bars. Values
%     of WIDTH > 1, produce overlapped bars.  The default value is WIDTH=0.8
%  
%     barh(...,'grouped') produces the default vertical grouped bar chart.
%     barh(...,'stacked') produces a vertical stacked bar chart.
%     barh(...,LINESPEC) uses the line color specified (one of 'rgbymckw').
%  
%     H = barh(...) returns a vector of handles to barseries objects.
%  
%     Use SHADING FACETED to put edges on the bars.  Use SHADING FLAT to
%     turn them off.
%  
%     See also hist, plot, barh, bar3, bar3h.
%     
%     In addition to those matlab properties, RISE adds further properties,
%     which allow to control for:
%     - the figure size (multiple plots): 'figsize', with default [3,3]
%     - the figure title (multiple plots): 'figtitle', with default ''
%     - the number of tick marks: 'nticks', with default 8
%     - the date format (see matlab's datestr) : 'date_format', with ''
%     - the log scale : 'logy', with default false
%     - the list of variables going to the secondary y axis : 'secondary_y'
%     - the flag for multiple plots: 'subplots' with default false

% bar('1994m7:1997m1',db(:,:,1),...
%     'figsize',[2,2],...
%     'figtitle','no title',...
%     'nticks',10,...
%     'legend',{'v1','v2'},...
%     'legend_loc','BO',...
%     'logy',true,...
%     'secondary_y',{'v1','v4'},...
%     'subplots',true,...
%     'linewidth',2,...
%     'date_format',17);
% xrotate(90)

plot_handle=utils.plot.myplot(@barh,varargin{:});

end