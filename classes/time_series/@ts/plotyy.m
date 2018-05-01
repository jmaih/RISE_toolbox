function varargout=plotyy(varargin)
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


%  plotyy Graphs with y tick labels on the left and right
%     plotyy(Y1,Y2) plots Y1 with y-axis labeling on the left and plots Y2
%     with y-axis labeling on the right. It uses the union of the dates of
%     the two time series to set the dates
%  
%     plotyy(xrange,Y1,Y2) uses the dates in xrange
%
%     plotyy(...,Y1,Y2,FUN) uses the plotting function FUN
%     instead of PLOT to produce each graph. FUN can be a
%     function handle or a string that is the name of a plotting
%     function, e.g. plot, semilogx, semilogy, loglog, stem,
%     etc. or any function that accepts the syntax H = FUN(X,Y).
%     For example
%        plotyy(...,Y1,Y2,@loglog)  % Function handle
%        plotyy(...,Y1,Y2,'loglog') % String
%  
%     plotyy(...,Y1,Y2,FUN1,FUN2) uses FUN1(...,Y1) to plot the data for
%     the left axes and FUN2(...,Y2) to plot the data for the right axes.
%  
%     [AX,H1,H2] = plotyy(...) returns the handles of the two axes created in
%     AX and the handles of the graphics objects from each plot in H1
%     and H2. AX(1) is the left axes and AX(2) is the right axes.
%
%     In addition to those matlab properties, RISE adds further properties,
%     which allow to control for. See parse_plot_args

% [ax,h1,h2]=plotyy('1994m7:1997m1',...
%     db(:,'v1',1),...
%     db(:,'v2',1),...
%     'nticks',10,...
%     'date_format',17)
%     xrotate(90)
    
[varargout{1:nargout}]=utils.plot.myplot(@plotyy,varargin{:});

end