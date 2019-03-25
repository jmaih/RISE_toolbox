% TEXT Add text descriptions to data points
%    TEXT(x,y,str) adds a text description to one or more data points in the
%    current axes using the text specified by str. To add text to one point,
%    specify x and y as scalars in data units. To add text to multiple
%    points, specify x and y as vectors with equal length.
% 
%    TEXT(x,y,z,str) positions the text in 3-D coordinates.
%    
%    TEXT(...,Name,Value) specifies text properties using one or more
%    Name,Value pair arguments. For example, 'FontSize',14 sets the font
%    size to 14 points. You can specify text properties with any of the
%    input argument combinations in the previous syntaxes. If you specify
%    the Position and String properties as Name,Value pairs, then you do not
%    need to specify the x, y, z, and str inputs.
% 
%    TEXT(container,...) creates the text in the axes, group, or transform
%    specified by container, instead of in the current axes.
% 
%    T = TEXT(...) returns one or more text objects. Use T to modify
%    properties of the text objects after they are created. For a list of
%    properties and descriptions, see Text Properties. You can specify an
%    output with any of the previous syntaxes.
% 
%    Execute GET(T), where T is a text object, to see a list of text object
%    properties and their current values.
%    Execute SET(T) to see a list of text object properties and legal
%    property values.
% 
%    See also XLABEL, YLABEL, ZLABEL, TITLE, GTEXT, LINE, PATCH.
%
%    Reference page in Doc Center
%       doc text
%
%