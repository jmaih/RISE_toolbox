% XLINE create a vertical line
%    XLINE(VALUES) creates constant lines at the specified values of x.
%    VALUE can be either a scalar or a vector. For example, xline(5) creates
%    a line at x = 5.  In 2-D views, the line is typically a vertical line.
%    In 3-D views, the line appears in the xy-plane in the middle of the
%    z-axis limits.
%    
%    XLINE(VALUES, LINESPEC) specifies either the line style, the color, or 
%    both. For example, ':' creates a dotted line, 'r' creates a red line, 
%    and 'r:' creates a dotted, red line.
% 
%    XLINE(VALUES, LINESPEC, LABELS) adds the specified labels to the lines.
%    When VALUES contains multiple elements, you can display different
%    labels for each line by specifying LABELS as either cell array of
%    character vectors or a string array with the same number of elements as
%    VALUES. Alternatively, specify LABELS as a character vector or a string
%    scalar to display the same label on all the lines.
% 
%    XLINE(...Name,Value) specifies ConstantLine properties using one or 
%    more name-value pair arguments. Specify name-value pairs after all 
%    other input arguments.
% 
%    XLINE(AX, ...) creates the line in the axes specified by ax instead of
%    in the current axes (gca). 
% 
%    XL = XLINE(...) returns the line. Use XL to modify the ConstantLine 
%    object after it is created.
% 
%    See also YLINE
%
%    Documentation for xline
%       doc xline
%
%    Other uses of xline
%
%       rise_dates.dates/xline
%