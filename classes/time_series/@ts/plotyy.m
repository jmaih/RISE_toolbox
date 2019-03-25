% PLOTYY Graphs with y tick labels on the left and right
%    PLOTYY(X1,Y1,X2,Y2) plots Y1 versus X1 with y-axis labeling
%    on the left and plots Y2 versus X2 with y-axis labeling on
%    the right.
% 
%    PLOTYY(X1,Y1,X2,Y2,FUN) uses the plotting function FUN
%    instead of PLOT to produce each graph. FUN can be a
%    function handle or a string that is the name of a plotting
%    function, e.g. plot, semilogx, semilogy, loglog, stem,
%    etc. or any function that accepts the syntax H = FUN(X,Y).
%    For example
%       PLOTYY(X1,Y1,X2,Y2,@loglog)  % Function handle
%       PLOTYY(X1,Y1,X2,Y2,'loglog') % String
% 
%    PLOTYY(X1,Y1,X2,Y2,FUN1,FUN2) uses FUN1(X1,Y1) to plot the data for
%    the left axes and FUN2(X2,Y2) to plot the data for the right axes.
% 
%    PLOTYY(AX,...) plots into AX as the main axes, instead of GCA.  If AX
%    is the vector of axes handles returned by a previous call to PLOTYY,
%    then the secondary axes will be ignored.
% 
%    [AX,H1,H2] = PLOTYY(...) returns the handles of the two axes created in
%    AX and the handles of the graphics objects from each plot in H1
%    and H2. AX(1) is the left axes and AX(2) is the right axes.
% 
%    See also PLOT, function_handle
%
%    Reference page in Doc Center
%       doc plotyy
%
%    Other functions named plotyy
%
%       ts/plotyy
%