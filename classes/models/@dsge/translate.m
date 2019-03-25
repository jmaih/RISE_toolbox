%--- help for polyshape/translate ---
%
%  TRANSLATE Translate a polyshape
% 
%  PG = TRANSLATE(pshape, V) translates a polyshape according to a 
%  two-element row vector V. The first element of V is the translation 
%  distance in the x direction, and the second element is the translation 
%  distance in the y direction. Positive values in V translate right and up, 
%  and negative values translate left and down. When pshape is an array of 
%  polyshapes, each element of pshape is translated according to V.
% 
%  PG = TRANSLATE(pshape, x, y) specifies the x and y translation distances 
%  as separate arguments.
% 
%  See also scale, rotate, polybuffer, polyshape
%
%    Reference page in Doc Center
%       doc polyshape/translate
%
%    Other functions named translate
%
%       dsge/translate
%