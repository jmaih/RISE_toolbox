% SLICE  Volumetric slice plot.
%    SLICE(X,Y,Z,V,Sx,Sy,Sz) draws slices along the x,y,z directions at
%    the points in the vectors Sx,Sy,Sz. The arrays X,Y,Z define the
%    coordinates for V and must be monotonic and 3-D plaid (as if
%    produced by MESHGRID).  The color at each point will be determined
%    by 3-D interpolation into the volume V.  V must be an M-by-N-by-P
%    volume array. 
% 
%    SLICE(X,Y,Z,V,XI,YI,ZI) draws slices through the volume V along the
%    surface defined by the arrays XI,YI,ZI.
% 
%    SLICE(V,Sx,Sy,Sz) or SLICE(V,XI,YI,ZI) assumes X=1:N, Y=1:M, Z=1:P. 
% 
%    SLICE(...,'method') specifies the interpolation method to use.
%    'method' can be 'linear', 'cubic', or 'nearest'.  'linear' is the
%    default (see INTERP3).
% 
%    SLICE(AX,...) plots into AX instead of GCA. 
% 
%    H = SLICE(...) returns a vector of handles to SURFACE objects.
% 
%    The axes CLim property is set to span the finite values of V.
% 
%    Example: To visualize the function x*exp(-x^2-y^2-z^2) over the
%    range -2 < x < 2, -2 < y < 2, -2 < z < 2, 
% 
%       [x,y,z] = meshgrid(-2:.2:2, -2:.25:2, -2:.16:2);
%       v = x .* exp(-x.^2 - y.^2 - z.^2);
%       slice(x,y,z,v,[-1.2 .8 2],2,[-2 -.2])
% 
%    See also MESHGRID, INTERP3.
%
%    Documentation for slice
%       doc slice
%
%