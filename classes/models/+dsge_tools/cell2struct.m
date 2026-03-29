% CELL2STRUCT Convert cell array to structure array.
%    S = CELL2STRUCT(C,FIELDS,DIM) converts the cell array C into
%    the structure S by folding the dimension DIM of C into fields of
%    S.  SIZE(C,DIM) must match the number of field names in FIELDS.
%    FIELDS can be a character array or a cell array of character vectors.
% 
%    Example:
%      c = {'tree',37.4,'birch'};
%      f = {'category','height','name'};
%      s = cell2struct(c,f,2);
% 
%    See also STRUCT2CELL, FIELDNAMES.
%
%    Documentation for cell2struct
%       doc cell2struct
%
%    Other uses of cell2struct
%
%       codistributed/cell2struct
%