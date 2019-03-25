% PARETO Pareto chart.
%    PARETO(Y,NAMES) produces a Pareto chart where the values in the
%    vector Y are drawn as bars in descending order.  Each bar will
%    be labeled with the associated name in the string matrix or
%    cell array NAMES.
% 
%    PARETO(Y,X) labels each element of Y with the values from X.
%    PARETO(Y) labels each element of Y with its index.
% 
%    PARETO(AX,...) plots into AX as the main axes, instead of GCA.
% 
%    [H,AX] = PARETO(...) returns a combination of patch and line object
%    handles in H and the handles to the two axes created in AX.
% 
%    See also HISTOGRAM, BAR.
%
%    Reference page in Doc Center
%       doc pareto
%
%