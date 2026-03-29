% PARETO Pareto chart.
%    PARETO(Y) produces a Pareto chart where the values in the vector Y
%    are drawn as bars in descending order.  Each bar will be labeled 
%    with with its element index in Y.
% 
%    PARETO(Y,NAMES) labels each bar with the associated text from the  
%    cell array or string vector NAMES.
% 
%    PARETO(Y,X) labels each bar of Y with the associated value from X.
% 
%    PARETO(...,THRESHOLD) Specify the amount of the cumulative
%    distribution displayed as a proportion between 0 and 1. The default
%    value for THRESHOLD is .95.
% 
%    PARETO(AX,...) plots into AX as the main axes, instead of GCA.
% 
%    H = PARETO(...) returns the primitive Line and Bar objects created.
% 
%    [H,AX] = PARETO(...) additionally returns the two axes objects 
%    created.
% 
%    See also HISTOGRAM, BAR.
%
%    Documentation for pareto
%       doc pareto
%
%