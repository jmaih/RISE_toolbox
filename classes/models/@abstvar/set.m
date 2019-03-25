% SET    Set object properties.
%    SET(H,'PropertyName',PropertyValue) sets the value of the specified
%    property for the graphics object with handle H.  H can be a vector of 
%    handles, in which case SET sets the properties' values for all objects
%    of H.
% 
%    SET(H,a) where a is a structure whose field names are object property
%    names, sets the properties named in each field name with the values
%    contained in the structure.
% 
%    SET(H,pn,pv) sets the named properties specified in the cell array of
%    strings pn to the corresponding values in the cell array pv for all
%    objects specified in H.  The cell array pn must be 1-by-N, but the cell
%    array pv can be M-by-N where M is equal to length(H) so that each
%    object will be updated with a different set of values for the list of 
%    property names contained in pn.
% 
%    SET(H,'PropertyName1',PropertyValue1,'PropertyName2',PropertyValue2,...)
%    sets multiple property values with a single statement.  Note that it
%    is permissible to use property/value string pairs, structures, and
%    property/value cell array pairs in the same call to SET.
% 
%    A = SET(H, 'PropertyName') 
%    SET(H,'PropertyName')
%    returns or displays the possible values for the specified property of 
%    the object with handle H.  The returned array is a cell array of 
%    possible value strings or an empty cell array if the property does not 
%    have a finite set of possible string values.
%    
%    A = SET(H) 
%    SET(H) 
%    returns or displays the names of the user-settable properties and 
%    their possible values for the object with handle H.  The return value 
%    is a structure whose field names are the user-settable property names 
%    of H, and whose values are cell arrays of possible property values or 
%    empty cell arrays.
% 
%    The default value for an object property can be set on any of an 
%    object's ancestors by setting the PropertyName formed by concatenating 
%    the string 'Default', the object type, and the property name.  For 
%    example, to set the default color of text objects to red in the current
%    figure window:
% 
%       set(gcf,'DefaultTextColor','red')
%    
%    Defaults can not be set on a descendant of the object, or on the
%    object itself - for example, a value for 'DefaultAxesColor' can not
%    be set on an axes or an axes child, but can be set on a figure or on
%    the root.
% 
%    Three strings have special meaning for PropertyValues:
%      'default' - use default value (from nearest ancestor)
%      'factory' - use factory default value
%      'remove'  - remove default value.
% 
%    See also GET, RESET, DELETE, GCF, GCA, FIGURE, AXES.
%
%    Reference page in Doc Center
%       doc set
%
%    Other functions named set
%
%       abstvar/set      editline/set      scribehandle/set
%       arrowline/set    editrect/set      scribehgobj/set
%       axischild/set    figobj/set        serial/set
%       axisobj/set      generic/set       splanar/set
%       axistext/set     gogetter/set      timer/set
%       COM/set          instrument/set    timeseries/set
%       dataset/set      RandStream/set    tscollection/set
%       dsge/set
%