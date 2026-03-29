%--- skipline.m not found. Showing help for spline instead. ---
%
% SPLINE Cubic spline data interpolation.
%    YQ = SPLINE(X,Y,XQ) performs cubic spline interpolation using the
%    values Y at sample points X to find interpolated values YQ at the query
%    points XQ.
%        - X must be a vector.
%        - If Y is a vector, Y(j) is the value at X(j).
%        - If Y is a matrix or n-D array, Y(:,...,:,j) is the value at X(j).
% 
%    SPLINE chooses slopes at X(j) such that YQ has a continuous second
%    derivative. Thus, SPLINE produces smooth results.
% 
%    Ordinarily, SPLINE uses not-a-knot conditions for the end slopes at
%    X(1) and X(end). However, if Y contains two more values than X has
%    entries, then the first and last value in Y are used as the end slopes.
% 
%    PP = SPLINE(X,Y) returns the piecewise polynomial form PP of the
%    interpolant. You can use PP as an input to PPVAL or UNMKPP.
% 
%    Comparison of SPLINE, PCHIP, and MAKIMA:
%        - All three are a form of piecewise cubic Hermite interpolation,
%          but each function computes the slopes of YQ at X(j) differently.
%        - SPLINE chooses slopes at X(j) such that the second derivative of
%          YQ is continuous. Therefore, SPLINE is smoother and more accurate
%          if the Y data represents values of a smooth function.
%        - PCHIP has no overshoots and less oscillation than SPLINE.
%        - MAKIMA has less oscillation than SPLINE but may have overshoots.
%        - PCHIP and MAKIMA are less expensive than SPLINE to set up PP.
%        - All three are equally expensive to evaluate.
%        - SPLINE and MAKIMA generalize to n-D grids. See INTERPN.
% 
%    Example: Compare SPLINE, PCHIP, and MAKIMA
% 
%        x = [1 2 3 4 5 5.5 7 8 9 9.5 10];
%        y = [0 0 0 0.5 0.4 1.2 1.2 0.1 0 0.3 0.6];
%        xq = 0.75:0.05:10.25;
%        yqs = spline(x,y,xq);
%        yqp = pchip(x,y,xq);
%        yqm = makima(x,y,xq);
% 
%        plot(x,y,'ko','LineWidth',2,'MarkerSize',10)
%        hold on
%        plot(xq,yqp,'LineWidth',4)
%        plot(xq,yqs,xq,yqm,'LineWidth',2)
%        legend('(x,y) data','pchip','spline','makima')
% 
%    Example: Interpolate a sine-like curve over a finer mesh
% 
%        x = 0:10;
%        y = sin(x);
%        xq = 0:.25:10;
%        yq = spline(x,y,xq);
%        figure
%        plot(x,y,'o',xq,yq)
% 
%    Example: Perform spline interpolation with prescribed end slopes.
%             Set the slopes to zero at the end points of the interpolant.
% 
%        x = -4:4;
%        y = [0 .15 1.12 2.36 2.36 1.46 .49 .06 0];
%        cs = spline(x,[0 y 0]);
%        xq = linspace(-4,4,101);
%        figure
%        plot(x,y,'o',xq,ppval(cs,xq));
% 
%    See also INTERP1, MAKIMA, PCHIP, PPVAL, MKPP, UNMKPP.
%
%    Documentation for spline
%       doc spline
%
%