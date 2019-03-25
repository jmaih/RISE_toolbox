% BAR Bar graph.
%    BAR(X,Y) draws the columns of the M-by-N matrix Y as M groups of N
%    vertical bars.  The vector X must not have duplicate values.
% 
%    BAR(Y) uses the default value of X=1:M.  For vector inputs, BAR(X,Y)
%    or BAR(Y) draws LENGTH(Y) bars.  The colors are set by the colormap.
% 
%    BAR(X,Y,WIDTH) or BAR(Y,WIDTH) specifies the width of the bars. Values
%    of WIDTH > 1, produce overlapped bars.  The default value is WIDTH=0.8
% 
%    BAR(...,'grouped') produces the default vertical grouped bar chart.
% 
%    BAR(...,'stacked') produces a vertical stacked bar chart.
% 
%    BAR(...,COLOR) uses the line color specified.  Specify the color as one of
%    these values: 'r', 'g', 'b', 'y', 'm', 'c', 'k', or 'w'.
% 
%    BAR(AX,...) plots into AX instead of GCA.
% 
%    H = BAR(...) returns a vector of handles to barseries objects.
% 
%    Examples: subplot(3,1,1), bar(rand(10,5),'stacked'), colormap(cool)
%              subplot(3,1,2), bar(0:.25:1,rand(5),1)
%              subplot(3,1,3), bar(rand(2,3),.75,'grouped')
% 
%    See also HISTOGRAM, PLOT, BARH, BAR3, BAR3H.
%
%    Reference page in Doc Center
%       doc bar
%
%    Other functions named bar
%
%       fints/bar    ts/bar
%