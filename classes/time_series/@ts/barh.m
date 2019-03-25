% BARH Horizontal bar graph.
%    BARH(X,Y) draws the columns of the M-by-N matrix Y as M groups of
%    N horizontal bars. The vector X must not have duplicate values.
% 
%    BARH(Y) uses the default value of X=1:M.  For vector inputs,
%    BARH(X,Y) or BARH(Y) draws LENGTH(Y) bars.  The colors are set by
%    the colormap.
% 
%    BARH(X,Y,WIDTH) or BARH(Y,WIDTH) specifies the width of the
%    bars. Values of WIDTH > 1, produce overlapped bars.  The
%    default value is WIDTH=0.8.
% 
%    BARH(...,'grouped') produces the default horizontal grouped bar chart.
%    BARH(...,'stacked') produces a horizontal stacked bar chart.
%    BARH(...,LINESPEC) uses the line color specified (one of 'rgbymckw').
% 
%    BARH(AX,...) plots into AX instead of GCA.
% 
%    H = BARH(...) returns a vector of handles to barseries objects.
% 
%    Use SHADING FACETED to put edges on the bars.  Use SHADING FLAT to
%    turn them off.
% 
%    Examples: subplot(3,1,1), barh(rand(10,5),'stacked'), colormap(cool)
%              subplot(3,1,2), barh(0:.25:1,rand(5),1)
%              subplot(3,1,3), barh(rand(2,3),.75,'grouped')
% 
%    See also PLOT, BAR, BAR3, BAR3H, HISTOGRAM.
% 
%    Copyright 1984-2017 The MathWorks, Inc.
%
%    Reference page in Doc Center
%       doc barh
%
%    Other functions named barh
%
%       fints/barh    ts/barh
%